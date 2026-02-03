"""Wrapper for running sequences through HuggingFace HIV-BERT models"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2024"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.0.2"

import os
import sys
import logging

# Configure logging to use Snakemake's log file
if "snakemake" in locals():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=snakemake.log[0] if hasattr(snakemake, 'log') and snakemake.log else None
    )
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

logger = logging.getLogger('hiv-bert-wrapper')

# Check if we need to use a custom environment
if "snakemake" in locals():
    custom_env = snakemake.params.get('custom_env', None)
    if custom_env:
        logger.info(f"Using custom environment: {custom_env}")
        # Add the custom environment to the Python path
        if os.path.exists(custom_env):
            if custom_env not in sys.path:
                sys.path.insert(0, custom_env)
            # Also add the site-packages directory
            site_packages = os.path.join(custom_env, 'lib', f'python{sys.version_info.major}.{sys.version_info.minor}', 'site-packages')
            if os.path.exists(site_packages) and site_packages not in sys.path:
                sys.path.insert(0, site_packages)
        else:
            logger.warning(f"Custom environment path does not exist: {custom_env}")

import pandas as pd
import numpy as np
from Bio import SeqIO, Seq
from Bio.Data import CodonTable
import pysam
import torch
from transformers import AutoModel, AutoTokenizer, AutoModelForSequenceClassification, pipeline
from contextlib import contextmanager
from pathlib import Path
import shutil
from typing import Optional, Tuple, Dict, Any
import yaml

if "snakemake" not in locals():
    import snakemake  # type: ignore

# Log CUDA information
logger.info(f"PyTorch version: {torch.__version__}")
logger.info(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    logger.info(f"CUDA version: {torch.version.cuda}")
    logger.info(f"GPU device: {torch.cuda.get_device_name(0)}")
    logger.info(f"GPU count: {torch.cuda.device_count()}")
else:
    logger.warning("CUDA is not available. Using CPU for inference.")

# HIV-BERT models
KNOWN_MODELS = {
    'damlab/hiv_bert': {'type': 'embedding'},
    'damlab/HIV_BERT': {'type': 'embedding'},
    'damlab/HIV_V3_bodysite': {'type': 'classification', 
                             'labels': ['periphery-tcell', 'periphery-monocyte', 'CNS', 
                                       'breast-milk', 'female-genitals', 'male-genitals', 
                                       'gastric', 'lung', 'organ']},
    'damlab/HIV_V3_coreceptor': {'type': 'classification',
                                'labels': ['CCR5', 'CXCR4']}
}

def translate_dna_to_aa(sequence, frame=0):
    """Translate DNA sequence to amino acids starting from specified frame"""
    if not sequence:
        return ""
    try:
        # Create a Seq object from the sequence
        dna_seq = Seq.Seq(sequence.upper())
        # Translate using the standard genetic code, starting from the specified frame
        aa_seq = str(dna_seq[frame:].translate(table=1, to_stop=False))
        return aa_seq
    except Exception as e:
        logger.warning(f"Translation error: {e} for sequence: {sequence[:20]}...")
        return ""

def format_seq_for_bert(sequence):
    """Format an amino acid sequence for BERT input"""
    # Replace rare amino acids with X
    for aa in 'UZOB':
        sequence = sequence.replace(aa, 'X')
    # Add spaces between amino acids
    sequence = ' '.join(sequence)
    return sequence

def get_model_path(model_name: str, model_directory: Optional[Path]) -> Optional[Path]:
    """Get path to cached model if it exists.
    
    Args:
        model_name: HuggingFace model name (e.g., 'damlab/hiv_bert')
        model_directory: Optional directory to cache models
        
    Returns:
        Path to cached model or None if not cached
    """
    if not model_directory:
        return None
        
    # Convert model name to directory structure (e.g., damlab/hiv_bert -> damlab_hiv_bert)
    safe_name = model_name.replace('/', '_')
    model_path = Path(model_directory) / safe_name
    
    # Check if model exists and has required files
    if model_path.exists():
        required_files = ['config.json', 'pytorch_model.bin', 'tokenizer.json']
        if all((model_path / file).exists() for file in required_files):
            return model_path
    
    return None

def save_model(model_name: str, model_path: Path, model: Any, tokenizer: Any) -> None:
    """Save model and tokenizer to specified directory.
    
    Args:
        model_name: HuggingFace model name
        model_path: Path to save directory
        model: Loaded model
        tokenizer: Loaded tokenizer
    """
    logger.info(f"Saving model to {model_path}")
    try:
        model.save_pretrained(model_path)
        tokenizer.save_pretrained(model_path)
    except Exception as e:
        logger.error(f"Failed to save model: {e}")
        # Clean up failed save attempt
        if model_path.exists():
            shutil.rmtree(model_path)
        raise

def load_model(model_name: str, model_directory: Optional[Path] = None) -> Tuple[Any, Any, Dict]:
    """Load model and tokenizer from HuggingFace or local cache.
    
    Args:
        model_name: HuggingFace model name
        model_directory: Optional directory to cache models
        
    Returns:
        Tuple of (model, tokenizer, model_info)
    """
    if model_name not in KNOWN_MODELS:
        logger.warning(f"Unknown model: {model_name}. Treating as general HuggingFace model.")
        model_info = {'type': 'unknown'}
    else:
        model_info = KNOWN_MODELS[model_name]
    
    # Check for cached model
    cached_path = None
    use_local = False
    if model_directory:
        cached_path = get_model_path(model_name, model_directory)
        if cached_path:
            logger.info(f"Loading model from cache: {cached_path}")
            use_local = True
            model_source = cached_path
        else:
            # Create directory for new model
            logger.info(f"Model not found in cache, downloading from HuggingFace: {model_name}")
            safe_name = model_name.replace('/', '_')
            cached_path = Path(model_directory) / safe_name
            cached_path.mkdir(parents=True, exist_ok=True)
            model_source = model_name
    else:
        logger.info(f"No cache directory specified, downloading from HuggingFace: {model_name}")
        model_source = model_name
    
    try:
        # Load tokenizer
        if use_local:
            tokenizer = AutoTokenizer.from_pretrained(model_source, local_files_only=True)
        else:
            tokenizer = AutoTokenizer.from_pretrained(model_source)
        
        # Load model based on type
        if model_info['type'] == 'embedding':
            if use_local:
                model = AutoModel.from_pretrained(model_source, local_files_only=True)
            else:
                model = AutoModel.from_pretrained(model_source)
        elif model_info['type'] == 'classification':
            if use_local:
                model = AutoModelForSequenceClassification.from_pretrained(model_source, local_files_only=True)
            else:
                model = AutoModelForSequenceClassification.from_pretrained(model_source)
        else:
            if use_local:
                model = AutoModel.from_pretrained(model_source, local_files_only=True)
            else:
                model = AutoModel.from_pretrained(model_source)
        
        # Save model to cache if it was downloaded
        if model_directory and not use_local:
            save_model(model_name, cached_path, model, tokenizer)
            
        return model, tokenizer, model_info
        
    except Exception as e:
        logger.error(f"Failed to load model: {e}")
        # Clean up failed cache attempt
        if cached_path and not use_local:
            if cached_path.exists():
                shutil.rmtree(cached_path)
        raise

def get_sequence_embedding(model, tokenizer, sequence, device='cpu'):
    """Get the embedding vector for a sequence using the base model"""
    formatted_seq = format_seq_for_bert(sequence)
    inputs = tokenizer(formatted_seq, return_tensors="pt").to(device)
    
    with torch.no_grad():
        outputs = model(**inputs)
    
    # Get the [CLS] token embedding (first token)
    embeddings = outputs.last_hidden_state[:, 0, :].cpu().numpy()[0]
    return embeddings

def get_sequence_classification(model, tokenizer, sequence, model_info, device='cpu'):
    """Get classification probabilities for a sequence"""
    formatted_seq = format_seq_for_bert(sequence)
    inputs = tokenizer(formatted_seq, return_tensors="pt").to(device)
    
    with torch.no_grad():
        outputs = model(**inputs)
    
    # Get probabilities
    logits = outputs.logits
    probabilities = torch.sigmoid(logits).cpu().numpy()[0]
    
    # Create a dictionary of label: probability
    result = {label: float(prob) for label, prob in zip(model_info['labels'], probabilities)}
    return result

def process_sequence(seq_id, aa_seq, model, tokenizer, model_info, device='cpu'):
    """Process a single sequence and return the result"""
    if model_info['type'] == 'embedding':
        embedding = get_sequence_embedding(model, tokenizer, aa_seq, device=device)
        result = {"id": seq_id}
        for j, value in enumerate(embedding):
            result[f"dim_{j}"] = value
        return result
    elif model_info['type'] == 'classification':
        probs = get_sequence_classification(model, tokenizer, aa_seq, model_info, device=device)
        result = {"id": seq_id, **probs}
        return result
    return None

def is_dna_sequence(sequence: str) -> bool:
    """Determine if a sequence is DNA based on composition.
    
    Args:
        sequence: Input sequence string
        
    Returns:
        True if sequence appears to be DNA
    """
    # Convert to uppercase for checking
    seq = sequence.upper()
    # Count DNA-specific bases
    dna_bases = sum(1 for base in seq if base in 'ATCGN')
    # If >80% of bases are DNA-specific, consider it DNA
    return (dna_bases / len(seq)) > 0.8 if seq else False

def process_sequence_string(sequence: str, frame: int = 0) -> str:
    """Process input sequence string to amino acids.
    
    Args:
        sequence: Input sequence (DNA or AA)
        frame: Reading frame if DNA sequence
        
    Returns:
        Amino acid sequence
    """
    if is_dna_sequence(sequence):
        return translate_dna_to_aa(sequence, frame)
    return sequence

def yield_sequences_from_fasta(file_path, params):
    """Generate sequences from a FASTA file"""
    logger.info(f"Processing FASTA file: {file_path}")
    for record in SeqIO.parse(file_path, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        
        # Process sequence (translate if DNA)
        aa_seq = process_sequence_string(sequence, params['translate_frame'])
        
        # Apply length filters
        if len(aa_seq) < params['min_length']:
            continue
        if params['max_length'] > 0 and len(aa_seq) > params['max_length']:
            aa_seq = aa_seq[:params['max_length']]
        
        yield seq_id, aa_seq

def yield_sequences_from_fastq(file_path, params):
    """Generate sequences from a FASTQ file"""
    logger.info(f"Processing FASTQ file: {file_path}")
    for record in SeqIO.parse(file_path, "fastq"):
        seq_id = record.id
        dna_seq = str(record.seq)
        
        # Translate DNA to amino acids
        aa_seq = translate_dna_to_aa(dna_seq, frame=params['translate_frame'])
        
        # Apply length filters
        if len(aa_seq) < params['min_length']:
            continue
        if params['max_length'] > 0 and len(aa_seq) > params['max_length']:
            aa_seq = aa_seq[:params['max_length']]
        
        yield seq_id, aa_seq

def yield_sequences_from_bam(file_path, params):
    """Generate sequences from a BAM/SAM file"""
    logger.info(f"Processing BAM/SAM file: {file_path}")
    with pysam.AlignmentFile(file_path, "rb") as bam:
        for read in bam:
            # Skip unmapped reads if required
            if params['mapped_only'] and read.is_unmapped:
                continue
                
            seq_id = read.query_name
            dna_seq = read.query_sequence
            
            if not dna_seq:
                continue
                
            # Translate DNA to amino acids
            aa_seq = translate_dna_to_aa(dna_seq, frame=params['translate_frame'])
            
            # Apply length filters
            if len(aa_seq) < params['min_length']:
                continue
            if params['max_length'] > 0 and len(aa_seq) > params['max_length']:
                aa_seq = aa_seq[:params['max_length']]
            
            yield seq_id, aa_seq

def process_sequences(sequence_generator, model, tokenizer, model_info, params):
    """Process sequences in batches from any generator"""
    results = []
    batch = []
    
    # Process sequences in batches
    for seq_id, aa_seq in sequence_generator:
        batch.append((seq_id, aa_seq))
        
        # Process batch when it reaches the target size
        if len(batch) >= params['batch_size']:
            for sid, seq in batch:
                result = process_sequence(sid, seq, model, tokenizer, model_info, params['device'])
                if result:
                    results.append(result)
            batch = []
    
    # Process any remaining sequences
    for seq_id, aa_seq in batch:
        result = process_sequence(seq_id, aa_seq, model, tokenizer, model_info, params['device'])
        if result:
            results.append(result)
    
    return results

def detect_file_format(file_path):
    """Detect the format of the input file"""
    if file_path.endswith('.fa') or file_path.endswith('.fasta'):
        return 'fasta'
    elif file_path.endswith('.fq') or file_path.endswith('.fastq'):
        return 'fastq'
    elif file_path.endswith('.bam'):
        return 'bam'
    elif file_path.endswith('.sam'):
        return 'sam'
    else:
        # Try to detect by content
        with open(file_path, 'rb') as f:
            header = f.read(4)
            if header.startswith(b'>'):
                return 'fasta'
            elif header.startswith(b'@'):
                with open(file_path, 'r') as text_file:
                    first_line = text_file.readline().strip()
                    second_line = text_file.readline().strip()
                    third_line = text_file.readline().strip()
                    if third_line.startswith('+'):
                        return 'fastq'
                    else:
                        return 'sam'
            elif header.startswith(b'\x1f\x8b'):  # gzip magic number
                return 'bam'
    
    raise ValueError(f"Could not determine file format for {file_path}")

def create_empty_output(model_info, output_file):
    """Create an empty output file with appropriate column headers"""
    if model_info['type'] == 'embedding':
        # Create empty DataFrame with id and dim_* columns
        df = pd.DataFrame(columns=['id'] + [f'dim_{i}' for i in range(768)])  # Typically 768 dimensions for BERT
    elif model_info['type'] == 'classification':
        # Create empty DataFrame with id and label columns
        df = pd.DataFrame(columns=['id'] + model_info['labels'])
    else:
        df = pd.DataFrame(columns=['id'])
        
    df.to_csv(output_file, index=False)
    logger.warning(f"No results to save. Created empty file: {output_file}")

def generate_metrics(results: list, model_info: dict) -> dict:
    """Generate metrics from processing results.
    
    Args:
        results: List of dictionaries containing processing results
        model_info: Dictionary containing model information
        
    Returns:
        Dictionary of metrics
    """
    if not results:
        return {
            'total_sequences': 0,
            'processed_sequences': 0,
            'mean_values': {},
            'high_confidence_counts': {}
        }
    
    df = pd.DataFrame(results)
    metrics = {
        'total_sequences': len(results),
        'processed_sequences': len(results)
    }
    
    # Calculate metrics based on model type
    if model_info['type'] == 'embedding':
        # For embedding model, calculate mean of each dimension
        embedding_cols = [col for col in df.columns if col.startswith('dim_')]
        metrics['mean_values'] = {
            col: float(df[col].mean()) 
            for col in embedding_cols
        }
        # No high confidence counts for embeddings
        metrics['high_confidence_counts'] = {}
        
    elif model_info['type'] == 'classification':
        # For classification model, calculate mean probability for each class
        prob_cols = [col for col in df.columns if col != 'id']
        metrics['mean_values'] = {
            col: float(df[col].mean()) 
            for col in prob_cols
        }
        # Count sequences with high confidence (>= 0.5) for each class
        metrics['high_confidence_counts'] = {
            col: int(df[df[col] >= 0.5].shape[0])
            for col in prob_cols
        }
    
    return metrics

def main():
    # Get input and output files
    input_file = str(snakemake.input[0])
    output_file = str(snakemake.output[0])
    metrics_file = snakemake.output.get('metrics', None)  # Get optional metrics file
    
    # Get parameters with defaults
    model_name = snakemake.params.get('model_name', 'damlab/hiv_bert')
    model_directory = snakemake.params.get('model_directory', None)
    if model_directory:
        model_directory = Path(model_directory)
    mapped_only = snakemake.params.get('mapped_only', False)
    translate_frame = snakemake.params.get('translate_frame', 0)
    min_length = snakemake.params.get('min_length', 10)
    max_length = snakemake.params.get('max_length', 256)
    batch_size = snakemake.params.get('batch_size', 32)
    sample_name = snakemake.params.get('sample_name', None)  # Add sample_name parameter
    
    # Determine device to use
    force_cpu = snakemake.params.get('force_cpu', False)
    if force_cpu:
        device = 'cpu'
        logger.info("Forcing CPU mode as requested")
    else:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        if device == 'cuda':
            logger.info(f"Using GPU: {torch.cuda.get_device_name(0)}")
        else:
            logger.warning("CUDA not available, falling back to CPU")
    
    params = {
        'mapped_only': mapped_only,
        'translate_frame': translate_frame,
        'min_length': min_length,
        'max_length': max_length,
        'batch_size': batch_size,
        'device': device
    }
    
    logger.info(f"Processing file: {input_file}")
    logger.info(f"Using model: {model_name}")
    logger.info(f"Parameters: {params}")
    
    try:
        # Load model and tokenizer
        model, tokenizer, model_info = load_model(model_name, model_directory)
        model.to(device)
        
        # Detect file format
        file_format = detect_file_format(input_file)
        logger.info(f"Detected file format: {file_format}")
        
        # Select sequence generator based on file format
        if file_format == 'fasta':
            sequence_generator = yield_sequences_from_fasta(input_file, params)
        elif file_format == 'fastq':
            sequence_generator = yield_sequences_from_fastq(input_file, params)
        elif file_format in ['bam', 'sam']:
            sequence_generator = yield_sequences_from_bam(input_file, params)
        else:
            raise ValueError(f"Unsupported file format: {file_format}")
        
        # Process sequences
        results = process_sequences(sequence_generator, model, tokenizer, model_info, params)
        
        # Save results to CSV
        if results:
            df = pd.DataFrame(results)
            df.to_csv(output_file, index=False)
            logger.info(f"Saved {len(results)} results to {output_file}")
        else:
            create_empty_output(model_info, output_file)
        
        # Generate and save metrics if requested
        if metrics_file:
            metrics = generate_metrics(results, model_info)
            # Add additional information to metrics
            metrics['model_name'] = model_name
            metrics['model_type'] = model_info['type']
            metrics['input_file'] = input_file
            metrics['parameters'] = params
            if sample_name:  # Add sample name if provided
                metrics['sample_name'] = sample_name
            
            with open(metrics_file, 'w') as f:
                f.write("# HIV-BERT Processing Metrics\n")
                yaml.dump(metrics, f, default_flow_style=False)
            logger.info(f"Saved metrics to {metrics_file}")
            
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        raise

if __name__ == "__main__":
    main() 
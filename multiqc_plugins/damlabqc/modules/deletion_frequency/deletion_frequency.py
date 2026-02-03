from collections import OrderedDict
import logging
import yaml

from multiqc import config
from multiqc.plots import bargraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.plots import violin

log = logging.getLogger('damlabqc.deletion_frequency')

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Deletion Frequency',
            anchor='deletion_frequency',
            href='',
            info='analyzes deletion frequencies from BAM files.'
        )

        # Initialize data storage
        self._deletion_freq_data = dict()
        self._regions = set()
        
        self._collect_log_files()
        self._add_general_stats()
        self._add_sections()
    
    def _collect_log_files(self):
        """Collect data from log files found by multiqc."""
        for f in self.find_log_files('deletion_frequency'):
            self._parse_deletion_freq_log(f)

        # Filter out ignored samples
        self._deletion_freq_data = self.ignore_samples(self._deletion_freq_data)

        # No samples found
        if len(self._deletion_freq_data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self._deletion_freq_data)} reports")
        
        # Write parsed report data
        self.write_data_file(self._deletion_freq_data, 'multiqc_deletion_frequency')
    
    def _parse_deletion_freq_log(self, f):
        """Parse single deletion frequency log file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            s_name = data.get('sample_name', f['fn'])
            s_name = self.clean_s_name(s_name, f)
            region = data.get('region_name', 'region')
            
            # Create nested structure if sample doesn't exist
            if s_name not in self._deletion_freq_data:
                self._deletion_freq_data[s_name] = {}
            
            # Store data under sample and region
            self._deletion_freq_data[s_name][region] = data
            self._regions.add(region)
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
            raise e
    
    def _add_general_stats(self):
        """Add to the general stats table at the top of the report."""
        general_stats_data = {}
        headers = OrderedDict()
        
        # Create headers for each region
        for region in sorted(self._regions):
            headers[f'partial_deletion_frequency_{region}'] = {
                'title': f'{region} Partial Del. Freq.',
                'description': f'Partial deletion frequency for {region}',
                'max': 1,
                'min': 0,
                'scale': 'RdYlBu-rev',
                'format': '{:,.2%}',
                'placement': 1000
            }
            headers[f'full_deletion_frequency_{region}'] = {
                'title': f'{region} Full Del. Freq.',
                'description': f'Full deletion frequency for {region}',
                'max': 1,
                'min': 0,
                'scale': 'RdYlBu-rev',
                'format': '{:,.2%}',
                'placement': 1100
            }
            headers[f'reads_covered_partial_deleted_{region}'] = {
                'title': f'{region} Partial Del. Reads',
                'description': f'Number of reads with partial deletion in {region}',
                'scale': 'PuRd',
                'format': '{:,.0f}',
                'hidden': True,
                'placement': 1200
            }
            headers[f'reads_covered_full_deleted_{region}'] = {
                'title': f'{region} Full Del. Reads',
                'description': f'Number of reads with full deletion in {region}',
                'scale': 'PuRd',
                'format': '{:,.0f}',
                'hidden': True,
                'placement': 1300
            }
            headers[f'reads_covering_required_{region}'] = {
                'title': f'{region} Reads Covering',
                'description': f'Number of reads covering the required region in {region}',
                'scale': 'YlGn',
                'format': '{:,.0f}',
                'hidden': True,
                'placement': 1400
            }
            headers[f'total_reads_{region}'] = {
                'title': f'{region} Total Reads',
                'description': f'Total number of reads in {region}',
                'scale': 'GnBu',
                'format': '{:,.0f}',
                'hidden': True,
                'placement': 1500
            }

        # Fill in the data
        for s_name, regions in self._deletion_freq_data.items():
            general_stats_data[s_name] = {}
            for region, data in regions.items():
                general_stats_data[s_name].update({
                    f'total_reads_{region}': data['total_reads'],
                    f'reads_covering_required_{region}': data['reads_covering_required'],
                    f'reads_covered_partial_deleted_{region}': data['reads_covered_partial_deleted'],
                    f'reads_covered_full_deleted_{region}': data['reads_covered_full_deleted'],
                    f'partial_deletion_frequency_{region}': data['partial_deletion_frequency'],
                    f'full_deletion_frequency_{region}': data['full_deletion_frequency']
                })
        
        self.general_stats_addcols(general_stats_data, headers)
    
    def _add_sections(self):
        """Add sections to MultiQC report."""
        self._add_deletion_violin()
        for region in sorted(self._regions):
            self._add_deletion_stats(region)
    
    def _add_deletion_violin(self):
        """Add violin plot of deletion frequencies for all regions."""
        # Create violin data for both partial and full deletions
        violin_data_partial = {}
        violin_data_full = {}
        
        # Collect deletion frequencies for each sample
        for s_name, regions in self._deletion_freq_data.items():
            violin_data_partial[s_name] = {}
            violin_data_full[s_name] = {}
            for region in sorted(self._regions):
                if region in regions:
                    violin_data_partial[s_name][region] = regions[region]['partial_deletion_frequency']
                    violin_data_full[s_name][region] = regions[region]['full_deletion_frequency']

        self.add_section(
            name='Partial Deletion Frequency Distribution',
            anchor='partial-deletion-frequency-violin',
            description='Distribution of partial deletion frequencies across samples for all regions',
            plot=violin.plot(violin_data_partial,
                            pconfig={
                                'id': 'partial_deletion_frequency_violin',
                                'title': 'Partial Deletion Frequency Distribution',
                                'ylab': 'Partial Deletion Frequency',
                                'ymin': 0,
                                'ymax': 1,
                                'violin_width': 0.8,
                                'violin_box': True
                            })
        )

        self.add_section(
            name='Full Deletion Frequency Distribution',
            anchor='full-deletion-frequency-violin',
            description='Distribution of full deletion frequencies across samples for all regions',
            plot=violin.plot(violin_data_full,
                            pconfig={
                                'id': 'full_deletion_frequency_violin',
                                'title': 'Full Deletion Frequency Distribution',
                                'ylab': 'Full Deletion Frequency',
                                'ymin': 0,
                                'ymax': 1,
                                'violin_width': 0.8,
                                'violin_box': True
                            })
        )
    
    def _add_deletion_stats(self, region):
        """Add deletion statistics plot for a specific region."""
        deletion_data = {}
        
        for s_name, regions in self._deletion_freq_data.items():
            if region in regions:
                d = regions[region]
                deletion_data[s_name] = {
                    'Total Reads': d['total_reads'],
                    'Reads Covering Required Region': d['reads_covering_required'],
                    'Reads with Partial Deletion': d['reads_covered_partial_deleted'],
                    'Reads with Full Deletion': d['reads_covered_full_deleted']
                }

        self.add_section(
            name=f'Deletion Statistics - {region}',
            anchor=f'deletion-frequency-stats-{region}',
            description=f'Distribution of reads with respect to partial and full deletions for {region}',
            plot=bargraph.plot(deletion_data,
                             pconfig={
                                 'id': f'deletion_frequency_stats_{region}',
                                 'title': f'Deletion Frequency: Read Statistics - {region}',
                                 'ylab': 'Number of Reads',
                                 'stacking': 'normal'
                             })
        ) 
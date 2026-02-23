from collections import OrderedDict
import logging
import yaml

from multiqc import config
from multiqc.plots import bargraph
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.deletion_block_detection')


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Deletion Block Detection',
            anchor='deletion_block_detection',
            href='',
            info='analyzes deletion blocks detected from BAM files.'
        )

        self._data = dict()
        
        self._collect_log_files()
        self._add_general_stats()
        self._add_sections()
    
    def _collect_log_files(self):
        """Collect data from log files found by multiqc."""
        for f in self.find_log_files('deletion_block_detection'):
            self._parse_log(f)

        self._data = self.ignore_samples(self._data)

        if len(self._data) == 0:
            raise ModuleNoSamplesFound
        
        log.info(f"Found {len(self._data)} reports")
        
        self.write_data_file(self._data, 'multiqc_deletion_block_detection')
    
    def _parse_log(self, f):
        """Parse single deletion block detection log file."""
        try:
            data = yaml.safe_load(f['f'])
            if data is None:
                log.debug(f"Could not parse file: {f['fn']}")
                return
                
            s_name = data.get('sample_name', f['fn'])
            s_name = self.clean_s_name(s_name, f)
            
            self._data[s_name] = data
            
        except Exception as e:
            log.warning(f"Error parsing file {f['fn']}: {e}")
    
    def _add_general_stats(self):
        """Add to the general stats table at the top of the report."""
        general_stats_data = {}
        
        for s_name, data in self._data.items():
            general_stats_data[s_name] = {
                'total_reads': data.get('total_reads', 0),
                'reads_with_deletions': data.get('reads_with_deletions', 0),
                'unique_deletion_count': data.get('unique_deletion_count', 0),
                'total_deletion_count': data.get('total_deletion_count', 0),
                'deletion_frequency': data.get('deletion_frequency', 0),
            }
        
        headers = OrderedDict()
        headers['deletion_frequency'] = {
            'title': 'Del. Freq.',
            'description': 'Fraction of reads containing at least one deletion',
            'max': 1,
            'min': 0,
            'scale': 'RdYlBu-rev',
            'format': '{:,.2%}',
            'placement': 1000
        }
        headers['reads_with_deletions'] = {
            'title': 'Reads w/ Del.',
            'description': 'Number of reads containing at least one deletion',
            'scale': 'PuRd',
            'format': '{:,.0f}',
            'placement': 1100
        }
        headers['unique_deletion_count'] = {
            'title': 'Unique Dels.',
            'description': 'Number of unique deletion blocks detected',
            'scale': 'OrRd',
            'format': '{:,.0f}',
            'placement': 1200
        }
        headers['total_deletion_count'] = {
            'title': 'Total Dels.',
            'description': 'Total number of deletions across all reads',
            'scale': 'YlOrRd',
            'format': '{:,.0f}',
            'hidden': True,
            'placement': 1300
        }
        headers['total_reads'] = {
            'title': 'Total Reads',
            'description': 'Total number of reads processed',
            'scale': 'GnBu',
            'format': '{:,.0f}',
            'hidden': True,
            'placement': 1400
        }
        
        self.general_stats_addcols(general_stats_data, headers)
    
    def _add_sections(self):
        """Add sections to MultiQC report."""
        self._add_deletion_bargraph()
        self._add_read_stats_bargraph()
    
    def _add_deletion_bargraph(self):
        """Add bar graph showing deletion counts per sample."""
        plot_data = {}
        
        for s_name, data in self._data.items():
            plot_data[s_name] = {
                'Unique Deletions': data.get('unique_deletion_count', 0),
                'Total Deletion Events': data.get('total_deletion_count', 0),
            }
        
        self.add_section(
            name='Deletion Counts',
            anchor='deletion-block-counts',
            description='Number of deletions detected per sample.',
            helptext='''
            This plot shows deletion statistics for each sample:
            
            * **Unique Deletions**: Number of distinct deletion blocks detected
            * **Total Deletion Events**: Total count of deletions across all reads (a deletion appearing in multiple reads is counted multiple times)
            ''',
            plot=bargraph.plot(
                plot_data,
                pconfig={
                    'id': 'deletion_block_counts',
                    'title': 'Deletion Block Detection: Deletion Counts',
                    'ylab': 'Count',
                    'cpswitch': True,
                    'cpswitch_counts_label': 'Counts',
                }
            )
        )
    
    def _add_read_stats_bargraph(self):
        """Add bar graph showing read statistics per sample."""
        plot_data = {}
        
        for s_name, data in self._data.items():
            total_reads = data.get('total_reads', 0)
            reads_with_dels = data.get('reads_with_deletions', 0)
            reads_without_dels = total_reads - reads_with_dels
            
            plot_data[s_name] = {
                'Reads with Deletions': reads_with_dels,
                'Reads without Deletions': reads_without_dels,
            }
        
        self.add_section(
            name='Read Statistics',
            anchor='deletion-block-read-stats',
            description='Distribution of reads with and without deletions.',
            helptext='''
            This plot shows the distribution of reads for each sample:
            
            * **Reads with Deletions**: Number of reads containing at least one deletion block
            * **Reads without Deletions**: Number of reads without any detected deletions
            
            Use the buttons above the plot to switch between absolute counts and percentages.
            ''',
            plot=bargraph.plot(
                plot_data,
                pconfig={
                    'id': 'deletion_block_read_stats',
                    'title': 'Deletion Block Detection: Read Statistics',
                    'ylab': 'Number of Reads',
                    'stacking': 'normal',
                    'cpswitch': True,
                    'cpswitch_counts_label': 'Number of Reads',
                    'cpswitch_percent_label': 'Percentage',
                }
            )
        )

from collections import OrderedDict
import logging
import yaml

from multiqc.plots import bargraph, heatmap, table
from multiqc.base_module import BaseMultiqcModule
from multiqc.base_module import ModuleNoSamplesFound

log = logging.getLogger('damlabqc.deletion_block_detection')


def _deletion_interval_key(d: dict) -> str:
    return f"{int(d['deletion_start'])}-{int(d['deletion_end'])}"


def _sort_interval_keys(keys: list) -> list:
    def key_tuple(k: str):
        a, _, b = k.partition('-')
        try:
            return (int(a), int(b))
        except ValueError:
            return (0, 0)

    return sorted(keys, key=key_tuple)


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

        def _has_targeted(d: dict) -> bool:
            tr = d.get('targeted_regions') or []
            tc = d.get('targeted_region_count')
            if tc is not None:
                return int(tc) > 0
            return len(tr) > 0

        any_targeted = any(_has_targeted(d) for d in self._data.values())

        for s_name, data in self._data.items():
            row = {
                'total_reads': data.get('total_reads', 0),
                'reads_with_deletions': data.get('reads_with_deletions', 0),
                'unique_deletion_count': data.get('unique_deletion_count', 0),
                'total_deletion_count': data.get('total_deletion_count', 0),
                'deletion_frequency': data.get('deletion_frequency', 0),
                'deletion_richness': data.get(
                    'deletion_richness', data.get('unique_deletion_count', 0)
                ),
                'deletion_shannon_entropy': data.get('deletion_shannon_entropy', 0.0),
            }
            if any_targeted:
                tr = data.get('targeted_regions') or []
                cov_sum = data.get('target_reads_covering_sum')
                if cov_sum is None and tr:
                    cov_sum = sum(int(x.get('reads_covering', 0)) for x in tr)
                elif cov_sum is None:
                    cov_sum = 0
                ovl_sum = data.get('target_reads_with_deletion_overlapping_sum')
                if ovl_sum is None and tr:
                    ovl_sum = sum(int(x.get('reads_with_deletion_overlapping', 0)) for x in tr)
                elif ovl_sum is None:
                    ovl_sum = 0
                row['target_reads_covering_sum'] = int(cov_sum)
                row['target_reads_with_deletion_overlapping_sum'] = int(ovl_sum)
                row['target_deletion_overlap_fraction'] = (
                    float(ovl_sum) / float(cov_sum) if cov_sum else 0.0
                )
            general_stats_data[s_name] = row

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
        headers['deletion_richness'] = {
            'title': 'Del. richness',
            'description': 'Distinct deletion block count (same as unique deletions; Shannon diversity input)',
            'scale': 'OrRd',
            'format': '{:,.0f}',
            'placement': 1250
        }
        headers['deletion_shannon_entropy'] = {
            'title': 'Del. Shannon H',
            'description': 'Shannon entropy (nats) of read counts across distinct deletion blocks',
            'scale': 'Viridis',
            'format': '{:,.3f}',
            'placement': 1275
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

        if any_targeted:
            headers['target_reads_covering_sum'] = {
                'title': 'Target cov. reads',
                'description': 'Sum over query regions of reads overlapping each region (same read in multiple regions counted multiple times)',
                'scale': 'Blues',
                'format': '{:,.0f}',
                'placement': 1450
            }
            headers['target_reads_with_deletion_overlapping_sum'] = {
                'title': 'Target del reads',
                'description': 'Sum over query regions of reads with a deletion overlapping each region',
                'scale': 'Purples',
                'format': '{:,.0f}',
                'placement': 1460
            }
            headers['target_deletion_overlap_fraction'] = {
                'title': 'Target del / cov.',
                'description': 'Ratio of summed target del-overlap reads to summed covering reads',
                'max': 1,
                'min': 0,
                'scale': 'RdYlBu-rev',
                'format': '{:,.2%}',
                'placement': 1470
            }

        self.general_stats_addcols(general_stats_data, headers)

    def _add_sections(self):
        """Add sections to MultiQC report."""
        self._add_deletion_bargraph()
        self._add_read_stats_bargraph()
        self._add_top_deletion_heatmap()
        self._add_top_deletion_table()

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

    def _union_top_deletion_keys(self) -> list:
        keys = set()
        for data in self._data.values():
            for d in data.get('top_deletions') or []:
                keys.add(_deletion_interval_key(d))
        return _sort_interval_keys(list(keys))

    def _read_count_for_key(self, data: dict, interval_key: str) -> int:
        ds, _, de = interval_key.partition('-')
        try:
            s0, s1 = int(ds), int(de)
        except ValueError:
            return 0
        for d in data.get('top_deletions') or []:
            if int(d['deletion_start']) == s0 and int(d['deletion_end']) == s1:
                return int(d['read_count'])
        return 0

    def _coverage_for_key(self, data: dict, interval_key: str):
        ds, _, de = interval_key.partition('-')
        try:
            s0, s1 = int(ds), int(de)
        except ValueError:
            return None
        for d in data.get('top_deletions') or []:
            if int(d['deletion_start']) == s0 and int(d['deletion_end']) == s1:
                return int(d['coverage_count'])
        return None

    def _add_top_deletion_heatmap(self):
        xcats = self._union_top_deletion_keys()
        if not xcats:
            return

        ycats = sorted(self._data.keys())
        heatmap_data = {}
        for s_name in ycats:
            heatmap_data[s_name] = {}
            d = self._data[s_name]
            for k in xcats:
                heatmap_data[s_name][k] = float(self._read_count_for_key(d, k))

        zmax = max(
            (v for row in heatmap_data.values() for v in row.values()),
            default=0.0,
        )
        if zmax <= 0:
            zmax = 1.0

        self.add_section(
            name='Top deletion blocks (heatmap)',
            anchor='deletion-block-top-heatmap',
            description='Read support for the union of each sample’s top 10 deletion intervals (reference coordinates). '
            'Columns are merged across samples; zeros mean the interval is not in that sample’s top ten.',
            helptext='''
            Intervals are labeled `start-end` on the reference. Values are **read_count** (how many reads carry that
            deletion block). Compare colors across samples to see shared vs sample-specific deletions.
            ''',
            plot=heatmap.plot(
                heatmap_data,
                xcats=xcats,
                ycats=ycats,
                pconfig={
                    'id': 'deletion_block_top_heatmap',
                    'title': 'Deletion Block Detection: Top deletions (read count)',
                    'xTitle': 'Deletion interval (start-end)',
                    'yTitle': 'Sample',
                    'min': 0,
                    'max': zmax,
                    'square': False,
                    'decimalPlaces': 0,
                    'legend': True,
                    'datalabels': True,
                }
            )
        )

    def _add_top_deletion_table(self):
        xcats = self._union_top_deletion_keys()
        if not xcats:
            return

        ycats = sorted(self._data.keys())
        table_data = OrderedDict()
        for s_name in ycats:
            d = self._data[s_name]
            total_reads = max(int(d.get('total_reads') or 0), 1)
            for k in xcats:
                rc = self._read_count_for_key(d, k)
                cov = self._coverage_for_key(d, k)
                rid = f'{s_name} | {k}'
                table_data[rid] = {
                    'sample': s_name,
                    'deletion_interval': k,
                    'read_count': rc,
                    'fraction_of_reads': float(rc) / float(total_reads) if total_reads else 0.0,
                    'coverage_count': int(cov) if cov is not None else 0,
                }

        headers = OrderedDict()
        headers['sample'] = {'title': 'Sample'}
        headers['deletion_interval'] = {'title': 'Interval (start-end)'}
        headers['read_count'] = {'title': 'Read count', 'format': '{:,.0f}'}
        headers['fraction_of_reads'] = {
            'title': 'Frac. total reads',
            'description': 'read_count / total_reads for this sample',
            'format': '{:,.4f}',
        }
        headers['coverage_count'] = {
            'title': 'Coverage',
            'description': 'Reads spanning the interval; 0 if interval not in that sample’s top-10 list',
            'format': '{:,.0f}',
        }

        self.add_section(
            name='Top deletion blocks (table)',
            anchor='deletion-block-top-table',
            description='Same matrix as the heatmap: one row per sample and per union interval.',
            helptext='''
            Sort and filter columns in the table UI. **Coverage** is only populated when that interval appears in
            that sample’s `top_deletions` list (top 10 by read count).
            ''',
            plot=table.plot(
                table_data,
                headers=headers,
                pconfig={
                    'id': 'deletion_block_top_table',
                    'title': 'Deletion Block Detection: Top deletions (table)',
                },
            )
        )

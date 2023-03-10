"""Script to generate haddock3 analysis html files from ../../test/golden_data directory.
"""
from pathlib import Path
from haddock.clis.cli_analyse import get_cluster_ranking
from haddock.libs.libplots import box_plot_handler, clt_table_handler, report_generator, scatter_plot_handler

golden_data = Path('../../tests/golden_data')
clt_file = golden_data / 'capri_clt_example.tsv'
ss_file = golden_data / 'capri_ss_example.tsv'
top_cluster = 10
format = None
scale = 1.0
step = str(golden_data)

html_files2keep = set(Path('.').glob('*.html'))

cluster_ranking = get_cluster_ranking(clt_file, top_cluster)
scatters = scatter_plot_handler(ss_file, cluster_ranking, format, scale)
boxes = box_plot_handler(ss_file, cluster_ranking, format, scale)
table = clt_table_handler(clt_file, ss_file)
report_generator(boxes, scatters, table, step)

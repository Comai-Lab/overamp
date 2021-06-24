# overamp

This is a python tool to remove overmaplified reads from NGS data. This verison takes a positon sorted BWA MEM mapped sam file and goes through checking for identical reads. It does this by pairs, so the read mush have the same starting location for both and non have differences within the read. In the case of nearly identical reads where one is a subsection of a larger, the larger will be kept. For run instructions please use ./overamp-7l-True-PE.py -h or ./overamp-7l-True-PE.py --help

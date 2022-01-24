#!/usr/bin/env python3

import argparse
import os
import re
import zipfile


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq-name', required=True)
    parser.add_argument('--fastqc-data-txt', required=True)
    parser.add_argument('--fastqc-report-html', required=True)
    parser.add_argument('--summary-txt', required=True)
    parser.add_argument('--output-directory', default=os.getcwd())
    args = parser.parse_args()

    #
    name = re.sub('.f(ast)?q.gz', '_fastqc', os.path.basename(args.fastq_name))
    output = os.path.join(args.output_directory, name + '.zip')

    with zipfile.ZipFile(output, 'w') as zf:
        zf.writestr(f'{name}/', '')
        zf.write(args.fastqc_data_txt, arcname=f'{name}/fastqc_data.txt')
        zf.write(args.fastqc_report_html, arcname=f'{name}/fastqc_report.html')
        zf.write(args.summary_txt, arcname=f'{name}/summary.txt')


if __name__ == '__main__':
    main()

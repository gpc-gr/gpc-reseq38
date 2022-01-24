#!/usr/bin/env python3

import argparse
import os

from bs4 import BeautifulSoup


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('source')
    parser.add_argument('output')
    parser.add_argument(
        '--assets-dir',
        default=os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(__file__))), 'assets')
    )
    args = parser.parse_args()

    #
    with open(args.source) as fin:
        soup = BeautifulSoup(fin.read(), 'html.parser')

    for tag in soup.find_all('link', rel='stylesheet'):
        with open(os.path.join(args.assets_dir, tag['href'].split('/')[-1])) as fin:
            new_tag = soup.new_tag('style')
            new_tag.append(fin.read())

        tag.replace_with(new_tag)

    for tag in soup.find_all('script'):
        if not tag.get('src'):
            continue

        with open(os.path.join(args.assets_dir, tag['src'].split('/')[-1])) as fin:
            new_tag = soup.new_tag('script')
            new_tag.append(fin.read())

        tag.replace_with(new_tag)

    #
    with open(args.output, 'w') as fout:
        fout.write(str(soup))


if __name__ == '__main__':
    main()

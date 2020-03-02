import os
import argparse
import yaml
import pprint
from easydict import EasyDict as edict

from download import download
from read_process import read_process
from de_analysis import de
from cancer import cancer

def parse_args():
    parser = argparse.ArgumentParser(description='eCLIP')
    parser.add_argument('--config', dest='config_file',
                        help='configuration filename',
                        default='configs.yml', type=str)
    return parser.parse_args()

def load_config(config_path):
    with open(config_path, 'r') as f:
        config = edict(yaml.load(f))
    return config

def main():

    print('ECLIP data processing pipeline.')

    ## load config file
    args = parse_args()
    if args.config_file is None:
        raise Exception('no configuration file')

    config = load_config(args.config_file)
    pprint.PrettyPrinter(indent=2).pprint(config)

    ## download data
    download(config)

    ## reads processing
    read_process(config)

    ## differential expression analysis
    de(config)

    ## cancer
    cancer(config)


if __name__ == '__main__':
    main()

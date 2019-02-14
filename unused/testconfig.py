import configparser, os, string


def tc():
    
    config = configparser.ConfigParser()
#    config.readfp(open(os.path.expanduser('~/cmb/cosmomc/data/CBI.dataset')))
    config.read('~/cmb/cosmomc/data/CBI.dataset')
    return config

def t2():
    return readparams(os.path.expanduser('~/cmb/cosmomc/data/CBI.dataset'))

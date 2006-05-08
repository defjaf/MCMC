import ClData

def main():
    files = ["data/MAXIMA.dataset", "data/B03_NA_21July05.newdat" ]

    datasets = [ ClData.ClDataSet(file) for file in files ]

    return datasets

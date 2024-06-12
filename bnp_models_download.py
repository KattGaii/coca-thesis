import benepar

if __name__ == '__main__':
    print("Downloading benepar_en3")
    benepar.download('benepar_en3')
    print("\nFinished downloading!")
    print("Downloading benepar_en3_large")
    benepar.download('benepar_en3_large')
    print("\nFinished downloading!")
    print("Finished all the downloads!")



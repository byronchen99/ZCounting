from pandas import DataFrame
import ROOT

def tree_to_df(tree, arrSize=5):
    #if tree has arrays with variable length, split into 'arrSize' new columns and fill empty values with 'NaN'
    df = DataFrame()
    for key in tree.dtype.names:
        if len(tree[key].shape) == 2:
            for i in range(tree[key].shape[1]):
                df[key+"_"+str(i)] = tree[key][:,i]
        else:
            df[key] = tree[key]
            # zero padding because variable length arrays in dataframe columns can not be stored properly
            for key in df.keys():
                if df[key].dtype =='O':
                    for i in range(arrSize):
                        df[key+"_"+str(i)] = df[key].apply(lambda x: x[i] if len(x) > i else float('NaN'))
                    df = df.drop([key],axis=1)
    return df

def to_RootTime(time, currentYear):
    # converts brilcalc time to root TDatime
    time =  time.split(" ")
    return ROOT.TDatime(currentYear, int(time[0].split("/")[0]),
                         int(time[0].split("/")[1]), int(time[1].split(":")[0]),
                         int(time[1].split(":")[1]), int(time[1].split(":")[2])).Convert() + 7200

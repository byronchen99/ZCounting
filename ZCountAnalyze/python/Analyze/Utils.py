from pandas import DataFrame


#write minitree into dataframe
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
import numpy as np
from hist import Hist


def templates_2d_to_1d(templates_2d, hPV): # templates_2d is a h5py file, hPV is a hist type object, and output is a dictionary of np arrays
    templates_1d = {}
    for template_type in templates_2d.keys():
        array_1d = np.matmul(np.array(templates_2d.get(template_type)), hPV.view())
    
        # hist = Hist.new.Regular(50, 66, 116, name = 'x').Double()
        # for i in range(len(array_1d)):
        #     hist[i] = array_1d[i]
        # templates_1d[template_type] = hist
        templates_1d[template_type] = array_1d

    return templates_1d
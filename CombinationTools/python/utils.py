import pdb


def simplify_uncertainties(inputFile, years=None, ratio=False, force_correlated=False, force_uncorrelated=False, prefix=""): 
    print("Load {0} and simplify uncertainties".format(inputFile))
    import sys, csv
    import numpy as np
    import uncertainties as unc

    # This function will combine the values for different years with uncertainties as specified in the input
    # text file. The file should have the format:
    #
    # Description,               Corr,2015,2016,2017,2018
    # Value,                     -,   4.21,40.99,49.79,67.86
    # Length scale,              C,   0.5, 0.8, 0.3, 0.2
    # Orbit drift,               C,   0.2, 0.1, 0.2, 0.1
    # ...
    # where the first line has the list of years
    # the second line has the list of values
    # and the all subsequent lines have, for each uncertainty, the correlation for that uncertainty and the value
    # (in %) of that uncertainty for each year.
    # The correlation should be 'C' for fully correlated uncertainties, 'U' for uncorrelated uncertainties, or for
    # partially correlated, P## where ## is the percent correlation (e.g. 'P70' for a 70% correlated uncertainty).

    correlations = {}
    uncertainties = {}
    with open(inputFile, "r") as csv_file:
        reader = csv.reader(csv_file, skipinitialspace=True)
        for i, row in enumerate(reader):
            if len(row) == 0 or row[0][0] == '#':
                continue
            systName=row[0]
            if i == 0:
                if years == None:
                    # if no years are given, take all years from the file 
                    years = row[2:]
                    select = [True for y in row[2:]]
                else:
                    new_years = []
                    for year in years:
                        if year not in row[2:]:
                            print("Error: year {0} not found in top line".format(year))
                            sys.exit(1)

                    # new array of years in correct order
                    new_years = []
                    for y in row[2:]:
                        if y in years:
                            new_years.append(y)
                    years = new_years

                    select = [y in years for y in row[2:]]
            elif i == 1:
                if systName not in ['Value', 'Luminosity']:
                    print("Error: expected first row to have vthe central alues")
                    sys.exit(1)
                values = np.array([float(x) for x in np.array(row[2:])[select]])
                if len(years) != len(values):
                    print("Error: number of values specified doesn't match number of years")
                    sys.exit(1)
            else:
                if row[1] != 'C' and row[1] != 'U' and row[1][0] != 'P':
                    print("Error: line {0}, correlation should be C, U, or P##".format(i))
                    sys.exit(1)
                correlations[systName] = row[1]
                uncertainties[systName] = np.array([float(x)/100 for x in np.array(row[2:])[select]])
                if len(years) != len(uncertainties[systName]):
                    print("Error: number of uncertainties for",systName,"doesn't match number of years")
                    sys.exit(1)


    # total_values = sum(values)
    # if ratio:
    #     total_values = values[1]/values[0]

    # total_uncertainty_sq = 0
    # for u in uncertainties:
    #     if (correlations[u] == 'C' or force_correlated) and not force_uncorrelated:
    #         # Correlated -- just add up individual uncertainties
    #         this_uncertainty = 0
    #         for i in range(len(years)):
    #             this_uncertainty += uncertainties[u][i]*values[i]
    #         if ratio:
    #             abs_uncertainty_0=uncertainties[u][0]*values[0]
    #             abs_uncertainty_1=uncertainties[u][1]*values[1]
    #             if abs_uncertainty_0>abs_uncertainty_1:
    #                 this_uncertainty = (values[1]/values[0])*(np.sqrt( (abs_uncertainty_1/values[1])**2 + (abs_uncertainty_0/values[0])**2 - 2*abs_uncertainty_1**2/(values[0]*values[1]) ))
    #             else:
    #                 this_uncertainty = (values[1]/values[0])*(np.sqrt( (abs_uncertainty_1/values[1])**2 + (abs_uncertainty_0/values[0])**2 - 2*abs_uncertainty_0**2/(values[0]*values[1]) ))
    #         total_uncertainty_sq += this_uncertainty**2
    #     elif correlations[u] == 'U' or force_uncorrelated:
    #         # Uncorrelated -- add in quadrature
    #         this_uncertainty_sq = 0
    #         for i in range(len(years)):
    #             this_uncertainty_sq += (uncertainties[u][i]*values[i])**2
    #         if ratio:
    #             abs_uncertainty_0=uncertainties[u][0]*values[0]
    #             abs_uncertainty_1=uncertainties[u][1]*values[1]
    #             this_uncertainty_sq = (values[1]/values[0])**2*( (abs_uncertainty_1/values[1])**2 + (abs_uncertainty_0/values[0])**2 )
    #         total_uncertainty_sq += this_uncertainty_sq
    #     elif correlations[u][0] == 'P':
    #         # Partially correlated
    #         frac_correlated = int(correlations[u][1:])/100.0
    #         this_uncertainty_corr = 0
    #         this_uncertainty_uncorr_sq = 0
    #         for i in range(len(years)):
    #             # Split up the SQUARED uncertainty into correlated and uncorrelated components.
    #             this_term_corr_sq = ((uncertainties[u][i]*values[i])**2)*frac_correlated
    #             this_term_uncorr_sq = ((uncertainties[u][i]*values[i])**2)*(1-frac_correlated)
    #             this_uncertainty_corr += np.sqrt(this_term_corr_sq)
    #             this_uncertainty_uncorr_sq += this_term_uncorr_sq
    #         total_uncertainty_sq += this_uncertainty_corr**2 + this_uncertainty_uncorr_sq

    # total_uncertainty = np.sqrt(total_uncertainty_sq)
    # if force_correlated:
    #     print("*** All uncertainties have been treated as correlated ***")
    # if force_uncorrelated:
    #     print("*** All uncertainties have been treated as uncorrelated ***")
    # if not ratio:
    #     print("Total values is %.2f +/- %.2f (uncertainty of %.2f%%)" % \
    #         (total_values, total_uncertainty, 100*total_uncertainty/total_values))
    # else:
    #     print("Ratio values is %.2f +/- %.3f (uncertainty of %.2f%%)" % \
    #         (total_values, total_uncertainty, 100*total_uncertainty/total_values))
                    
    # # Next make the final table for use by other people. The first step is to go through and see if any
    # # uncertainties are treated as correlated but only are nonzero for one year. If so, we can treat them as
    # # uncorrelated rather than have to put them as a separate bucket.
    # for u in uncertainties:
    #     n_nonzero = 0
    #     for x in uncertainties[u]:
    #         if x > 0:
    #             n_nonzero += 1
    #     if n_nonzero == 1 and correlations[u] == 'C':
    #         correlations[u] = 'U'

    # # Now, for each uncertainty, print it out as is if correlated, but add it to the total uncorrelated bin if
    # # not.
    # total_uncorrelated = [0]*len(years)
    # for u in sorted(uncertainties):
    #     if correlations[u] == 'U':
    #         for i in range(len(uncertainties[u])):
    #             total_uncorrelated[i] += (100*uncertainties[u][i])**2
    #     elif correlations[u][0] == 'P':
    #         correlated_uncertainty = []
    #         frac_correlated = int(correlations[u][1:])/100.0
    #         for i in range(len(years)):
    #             # Split the uncertainty into correlated and uncorrelated parts.
    #             this_term_corr_sq = ((100*uncertainties[u][i])**2)*frac_correlated
    #             this_term_uncorr_sq = ((100*uncertainties[u][i])**2)*(1-frac_correlated)
    #             # Dump the uncorrelated part into the uncorrelated bucket and keep track of the correlated part.
    #             total_uncorrelated[i] += this_term_uncorr_sq
    #             correlated_uncertainty.append("%.1f" % (np.sqrt(this_term_corr_sq)))
    # # Finally, print out the uncorrelateds.
    # for i in range(len(years)):
    #     output_array = ['0.0']*len(years)
    #     output_array[i] = "%.1f" % (np.sqrt(total_uncorrelated[i]))

    # build up covariance matrix

    # diagnoal matrix for uncorrelated uncertainties
    uncorr = np.zeros((len(years),len(years)))
    for i in range(len(years)):
        uncorr[i,i] = 1.0

    #covariance matrix to be filled
    covariance = np.zeros((len(years),len(years)))

    result = {}
    for key, uncert in uncertainties.items():
        #check which years are correlated
        if correlations[key] == 'U':
            for i, y in enumerate(years):
                uname = prefix+y
                if uncert[i] == 0: # skip if uncertainty is 0
                    continue
                if uname not in result:
                    result[uname] = {y: 0}
                result[uname][y] += uncert[i]**2
            
            covariance += uncorr * (uncert*values)**2

        elif correlations[key] == 'C':
            thisSet=[]
            for i in range(len(years)):
                if uncert[i] > 0:
                    thisSet.append(years[i])
                
            if len(thisSet) == 0:
                continue

            setName = prefix+"".join(thisSet)
            if setName not in result.keys():
                result[setName] = {y:0 for y in thisSet}

            for i, y in enumerate(years):
                if uncert[i] > 0:
                    result[setName][y] += uncert[i]**2

            corr = np.ones((len(years),len(years)))
            for i in range(len(values)):
                for j in range(len(values)):
                    corr[i,j] = uncert[i]  * values[i] * uncert[j] * values[j]

            covariance += corr

    values = unc.correlated_values(values, covariance)

    for unc_name, unc_values in result.items():
        for unc_year, unc_value in unc_values.items():
            result[unc_name][unc_year] = np.sqrt(unc_value)*100

    values = {y:v for y,v in zip(years, values)}

    return result, values


import numpy as np
import ROOT
import matplotlib.pyplot as plt
import uncertainties as unc

lPHYS = np.array([36.330, 41.480, 59.830])
lPHYS_Unc = np.array([0.424, 0.962, 1.501])

cPHYS =  np.array([
    [1.00,0.20,0.41],
    [0.20,1.00,0.34],
    [0.41,0.34,1.00]
])

uPHYS = unc.correlated_values_norm(zip(lPHYS, lPHYS_Unc), cPHYS)

# --- from combination of high pileup Z counts
hasZ = [0,1,1,0]
lZ = np.array([41.480,59.830])
lZ_Unc = np.array([0.637,1.486])

cZ =  np.array([
    [1.00,1.00],
    [1.00,1.00]
])

uZ = unc.correlated_values_norm(zip(lZ, lZ_Unc), cZ)


lComb = np.array([36.330,41.480,59.830])
lComb_Unc = np.array([0.405,0.536,1.060])

cComb =  np.array([
    [1.00,0.81,0.56],
    [0.81,1.00,0.94],
    [0.56,0.94,1.00]
])

uComb = unc.correlated_values_norm(zip(lComb, lComb_Unc), cComb)

# # --- from combination of 2017H lumi with high pileup Z counts
# hasZ = [1,1,1,1]
# lZ = np.array([36.330,41.480,59.830])
# lZ_Unc = np.array([0.742,0.847,1.222])
#
# cZ =  np.array([
#     [1.00,1.00,1.00],
#     [1.00,1.00,1.00],
#     [1.00,1.00,1.00]
# ])
#
# uZ = unc.correlated_values_norm(zip(lZ, lZ_Unc), cZ)
#
#
# lComb = np.array([36.330,41.480,59.830])
# lComb_Unc = np.array([0.397,0.464,0.672])
#
# cComb =  np.array([
#     [1.00,0.97,0.97],
#     [0.97,1.00,0.97],
#     [0.97,0.97,1.00]
# ])
#
# uComb = unc.correlated_values_norm(zip(lComb, lComb_Unc), cComb)



textsize=15


fig, ax = plt.subplots(1, 4, sharey=True, figsize=(10, 3))
fig.subplots_adjust(left=0.1, right=0.97, top=0.99, bottom=0.25, wspace=0.0)

xZOfffset = 0
for i in range(3):

    if hasZ[i] == 0:
        xZOfffset += 1
        xx = np.array([uPHYS[i].n, uComb[i].n])
        xxErr = np.array([uPHYS[i].s, uComb[i].s])
        yy = np.array([3., 1.])
        yyErr = np.array([0., 0.])
    else:
        yy = np.array([3., 2., 1.])
        yyErr = np.array([0., 0., 0.])
        xx = np.array([uPHYS[i].n, uZ[i-xZOfffset].n, uComb[i].n])
        xxErr = np.array([uPHYS[i].s, uZ[i-xZOfffset].s, uComb[i].s])

    ax[i].bar(uComb[i].n, 2.75, uComb[i].s*2, 0.5, color='blue', alpha=0.5)

    ax[i].errorbar(xx, yy, xerr=xxErr, yerr=yyErr, fmt='.k', ms=8, capsize=3)

if hasZ[3] == 0:
    yy = np.array([3., 1.])
    yyErr = np.array([0., 0.])
    xx = np.array([sum(uPHYS).n, sum(uComb).n])
    xxErr = np.array([sum(uPHYS).s, sum(uComb).s])
else:
    yy = np.array([3., 2., 1.])
    yyErr = np.array([0., 0., 0.])
    xx = np.array([sum(uPHYS).n, sum(uZ).n, sum(uComb).n])
    xxErr = np.array([sum(uPHYS).s, sum(uZ).s, sum(uComb).s])

ax[3].bar(sum(uComb).n, 2.75, sum(uComb).s*2, 0.5, color='blue', alpha=0.5)
ax[3].errorbar(xx, yy, xerr=xxErr, yerr=yyErr, fmt='.k', ms=8, capsize=3)



ax[0].text(0.1, 0.85, '2016', color='black', transform=ax[0].transAxes, fontsize=textsize)
ax[1].text(0.1, 0.85, '2017', color='black', transform=ax[1].transAxes, fontsize=textsize)
ax[2].text(0.1, 0.85, '2018', color='black', transform=ax[2].transAxes, fontsize=textsize)
ax[3].text(0.1, 0.85, 'RunII', color='black', transform=ax[3].transAxes, fontsize=textsize)
ax[3].text(0.6, 0.85, 'CMS', color='black', transform=ax[3].transAxes, fontsize=textsize*1.5, weight='bold')


ax[0].set_ylim(0.5, 4.)
labels = ['Comb.', 'Z Boson', 'PHYSICS']
ax[0].set_yticks([1., 2., 3.])
ax[0].set_yticklabels(labels, size=textsize)

for i in range(4):
    ax[i].yaxis.set_ticks_position('none')
    ax[i].tick_params(axis='x', labelsize=textsize)

ax[3].set_xlabel('Luminosity [fb$^{-1}$]', fontsize = textsize)

plt.savefig('LumiCombination.png')
plt.savefig('LumiCombination.pdf')

plt.close()

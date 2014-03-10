from pylab import *

pctidentities = []
nconflicts = []
with open('aln.txt', 'r') as infile:
    for l, line in enumerate(infile.readlines()):
        if l % 3 == 1:
            words = line.split()
            pctidentities.append( float(words[1]) )
            nconflicts.append( int(words[2].split('/')[0]) )

pctidentities = array(pctidentities)
nconflicts = array(nconflicts)
pctconflicts = 100. - pctidentities

nbins = 20
#hist(log10(nconflicts+1), nbins)
#hist(log10(pctconflicts+1), nbins)
hist(pctconflicts, nbins, normed=True, log=True)

xlabel('% conflicts')
ylabel('log_10 density')
ax = gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
#xlim(0,100)
#ylim(0,100)

# summary statistics
print 'Number of plasmids with 0 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts == 0), len(nconflicts), (sum(nconflicts == 0)*100./len(nconflicts)))
print 'Number of plasmids with 1 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts == 1), len(nconflicts), (sum(nconflicts == 1)*100./len(nconflicts)))
print 'Number of plasmids with 2 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts == 2), len(nconflicts), (sum(nconflicts == 2)*100./len(nconflicts)))
print 'Number of plasmids with > 0 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 0), len(nconflicts), (sum(nconflicts > 0)*100./len(nconflicts)))
print 'Number of plasmids with > 1 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 1), len(nconflicts), (sum(nconflicts > 1)*100./len(nconflicts)))
print 'Number of plasmids with > 2 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 2), len(nconflicts), (sum(nconflicts > 2)*100./len(nconflicts)))
print 'Number of plasmids with > 3 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 3), len(nconflicts), (sum(nconflicts > 3)*100./len(nconflicts)))
print 'Number of plasmids with > 10%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 10.), len(nconflicts), (sum(pctconflicts > 10.)*100./len(nconflicts)))
print 'Number of plasmids with > 25%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 25.), len(nconflicts), (sum(pctconflicts > 25.)*100./len(nconflicts)))
print 'Number of plasmids with > 50%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 50.), len(nconflicts), (sum(pctconflicts > 50.)*100./len(nconflicts)))
print 'Number of plasmids with > 75%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 75.), len(nconflicts), (sum(pctconflicts > 75.)*100./len(nconflicts)))

savefig('hist.png')


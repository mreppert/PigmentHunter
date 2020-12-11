## run through pymol, eg.:
## pymol -qc mutate.py -- 1god A/94/ ASN

from pymol import cmd
import sys

fname, selection, mutant = sys.argv[-3:]
cmd.wizard("mutagenesis")
cmd.load(fname)
cmd.refresh_wizard()
cmd.get_wizard().do_select(selection)
cmd.get_wizard().set_mode(mutant)

rotid = -1
best = 1e6
scores = cmd.get_wizard().bump_scores
for n in range(0, len(scores)):
    score = scores[n]
    print(score)
    if score < best:
        best = score
        rotid = n + 1

if len(scores)>0:
    print("****************************************************")
    print("Best score: "+str(best)+" for rotamer #" + str(rotid))
    cmd.frame(rotid)

cmd.get_wizard().apply()
cmd.set_wizard()
cmd.save(fname+".mut")
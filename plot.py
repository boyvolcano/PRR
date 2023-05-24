from sympy import symbols, exp, log
from matplotlib import pyplot as plt
import os

a = symbols('a')
n = symbols('(n*)')
pwd = os.getcwd()

evs = [
	('QuickSelect', exp(2 * a - a * log(a)), exp(1.15-0.28*a)),
	('QuickSort', exp((4 - a) * log(n)), exp(0.5-0.5*a)),
	('L1Diameter', exp(a - a * log(a)), exp(1.39-0.69*a)),
   	('L2Diameter', exp(a - a * log(a)), exp(1.39-0.69*a)),
	('RandSearch', exp((2 * a - a * log(a)) * log(n)), exp(-0.29*a*log(n))),
	('Channel', exp(0.5 * (2-a) * n), exp(1-0.37*a)),
	('Rdwalk', exp(0.25 * (1-a) * n), exp(0.6-0.41*a)),
	('Rdadder', exp(0.25 * (1-a) * n), exp(0)),
	('MC1', exp((a - a * log(a)) * log(n)), exp(-0.69*a*log(n))),
	('MC2', exp((a - a * log(a)) * log(n)), exp(-0.69*a*log(n))),
	('MC3', exp(2 * a - a * log(a)), exp(1.15-0.28*a)),
	('MC4', exp(a - a * log(a)), exp(0))]

sub_list = [{a:10,n:13}, {a:11, n:15}, {a:12, n:17}]
psub_list = [{a:i*0.25, n:17} for i in range(40,61)]
ulist = [i * 0.25 for i in range(40,61)]

def print_map(mps):
	return '$\\alpha={0};n^*={1}$'.format(mps[a], mps[n])

def transform(fnum):
	fnum = str(fnum)
	a = fnum[0:fnum.find('E')]
	b = fnum[fnum.find('E')+1:]
	if(b=='0' or b=='+0'):
		return '$' + str(a) + '$'
	if(b[0]=='+'):
		b = b[1:]
	return str('${0}*10^{{{1}}}$'.format(a,b))

print('Benchmark,Concrete choice,Our bound,Karp Bound,Ratio')

for (name, our_bound, karp_bound) in evs:
	for csub in sub_list:
		our_concrete = our_bound.evalf(subs=csub)
		karp_concrete = karp_bound.evalf(subs=csub)
		ratio = karp_concrete / our_concrete

		our_concrete = transform(format(our_concrete, '.2E'))
		karp_concrete = transform(format(karp_concrete, '.2E'))
		ratio = transform(format(ratio, '.2E'))

		if(karp_concrete == '$1.00$'):
			karp_concrete = 'Not applicable'
			ratio = '-'

		print(str(name) + ',' + print_map(csub) + ',' + str(our_concrete) + ',' + str(karp_concrete) + ',' + str(ratio))

idx = 1
for (name, our_boudn, karp_bound) in evs:
	print('\\begin{{figure}}\n\t\\centering\n\t\\includegraphics[width=0.7\\textwidth]{{figs/{0}.png}}\n\t\\caption{{The Plot for {0}}}\n\t\\label{{fig:b{2}}}\n\\end{{figure}}'.format(name,pwd,idx))
	idx = idx + 1


for (name, our_boudn, karp_bound) in evs:
	print('<img src=\"{1}/fig/{0}.png\" alt=\"{0}\" style=\"zoom:80%\" />\n<center>The plot for {0}</center>'.format(name,pwd))

for (name, our_bound, karp_bound) in evs:
	our_result = [(our_bound).evalf(subs=i) for i in psub_list]
	karp_result = [(karp_bound).evalf(subs=i) for i in psub_list]
	plt.plot(ulist, our_result, label = 'Our bound')
	if(not (name in ['Rdadder', 'MC4'])):
		plt.plot(ulist, karp_result, label = 'Karp\'s bound')
	plt.xlabel('The choice of α')
	plt.ylabel('The tail bound')
	plt.yscale('log')
	plt.legend()
	plt.savefig('fig/{0}.png'.format(name))
	plt.clf()


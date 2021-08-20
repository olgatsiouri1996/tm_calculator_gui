# python3
from gooey import *
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
# imput parameters
@Gooey(required_cols=2, program_name='Tm calculator', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-forward", "--forward", required=True, type=str, help="forward primer sequence")
	ap.add_argument("-reverse", "--reverse", required=True, type=str, help="forward primer sequence")
	ap.add_argument("-progam", "--program", required=False, type=int, default=1, help="program to choose(1.classic Tm calculation, 2.Calculation based on nearest neighbor thermodynamics). Default is 1")
	ap.add_argument("-Na", "--Na", required=False, type=float, default=50.0, help="Na+ concentration(float number)")
	ap.add_argument("-Mg", "--Mg", required=False, type=float, default=1.5, help="Mg2+ concentration(float number)")
	ap.add_argument("-dNTPs", "--dNTPs", required=False, type=float, default=0.6, help="dNTPs concentration(float number)")
	args = vars(ap.parse_args())
# create the function to calculate the Tm:
	def tm(seq):
		return 2 * seq.count("A") + 2 * seq.count("C") + 4 * seq.count("G") + 4 * seq.count("C")
# main
# import forward primer
	frwd = Seq(args['forward'])
# import reverse primer
	rever = Seq(args['reverse'])
# choose program
	program = args['program']
	if program == 1:
# calculate Tm of forward primer
		print("The Tm of the forward primer is:",tm(frwd),sep=" ")
# calculate Tm of reverse primer
		print("The Tm of the reverse primer is:",tm(rever),sep=" ")
	elif program == 2:
# calculate Tm of forward primer
		print("The Tm of the forward primer is:",'%0.2f' % mt.Tm_NN(frwd, Na=50, Mg=1.5, dNTPs=0.6),sep=" ")
# calculate Tm of reverse primer
		print("The Tm of the reverse primer is:",'%0.2f' % mt.Tm_NN(rever, Na=50, Mg=1.5, dNTPs=0.6),sep=" ")

if __name__ == '__main__':
	main()

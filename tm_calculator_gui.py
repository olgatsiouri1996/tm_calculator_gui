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
	ap.add_argument("-Na", "--Na", required=False, type=float, default=50.0, help="Na+ concentration")
	ap.add_argument("-Mg", "--Mg", required=False, type=float, default=1.5, help="Mg2+ concentration")
	ap.add_argument("-dNTPs", "--dNTPs", required=False, type=float, default=0.6, help="dNTPs concentration")
	args = vars(ap.parse_args())
# import forward primer
	frwd = Seq(args['forward'])
# import reverse primer
	rever = Seq(args['reverse'])
# calculate Tm of forward primer
	print("The Tm of the forward primer is:",round(mt.Tm_NN(frwd, Na=args['Na'], Mg=args['Mg'], dNTPs=args['dNTPs']), 2),sep=" ")
# calculate Tm of reverse primer
	print("The Tm of the reverse primer is:",round(mt.Tm_NN(rever, Na=args['Na'], Mg=args['Mg'], dNTPs=args['dNTPs']), 2),sep=" ")

if __name__ == '__main__':
	main()

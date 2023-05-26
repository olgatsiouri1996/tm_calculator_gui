# python3
from gooey import *
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
# imput parameters
@Gooey(required_cols=2, program_name='Tm calculator', header_bg_color= '#DCDCDC', terminal_font_color= '#000000', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-forward", "--Forward primer", required=True, type=str, help="forward primer sequence")
	ap.add_argument("-reverse", "--Reverse primer", required=True, type=str, help="forward primer sequence")
	ap.add_argument("-dmso", "--DSMO?", required=False, type=str, default="no", choices=["yes","no"], help="correct Tm calculation based on DMSO concertrations")
	ap.add_argument("-per", "--%DSMO", required=False, type=float, default=2.0, help="%DMSO concertration")
	ap.add_argument("-decrease", "--decrease of DMSO%", required=False, type=float, default=0.65, help="How much should Tm decreases per percent DMSO. Default=0.65 (von Ahsen et al. 2001)")
	ap.add_argument("-Na", "--Na", required=False, type=float, default=50.0, help="Na+ concentration")
	ap.add_argument("-Mg", "--Mg", required=False, type=float, default=1.5, help="Mg2+ concentration")
	ap.add_argument("-dNTPs", "--dNTPs", required=False, type=float, default=0.6, help="dNTPs concentration")
	args = vars(ap.parse_args())
# import forward primer
	frwd = Seq(args['Forward primer'])
# import reverse primer
	rever = Seq(args['Reverse primer'])
# calculate Tm of forward primer
	tmf = round(mt.Tm_NN(frwd, Na=args['Na'], Mg=args['Mg'], dNTPs=args['dNTPs']), 2)
# calculate Tm of reverse primer
	tmr = round(mt.Tm_NN(rever, Na=args['Na'], Mg=args['Mg'], dNTPs=args['dNTPs']), 2)
# if DMSO is added fix Tm calculation
	if args['DMSO?'] == "no":
		print(f'The Tm of the forward primer is:{tmf}',sep=" ")
		print(f'The Tm of the reverse primer is:{tmr}',sep=" ")
	else:
		fixedtmf = round(mt.chem_correction(tmf, DMSO=args['%DMSO'],DMSOfactor=args['decrease of DMSO%']) ,2)
		print(f'The Tm of the forward primer is:{fixedtmf}',sep=" ")
		fixedtmr = round(mt.chem_correction(tmr, DMSO=args['%DMSO'],DMSOfactor=args['decrease of DMSO%']) ,2)
		print(f'The Tm of the forward primer is:{fixedtmr}',sep=" ")

if __name__ == '__main__':
	main()

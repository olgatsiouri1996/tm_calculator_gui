# python3
from gooey import *
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
# imput parameters
@Gooey(required_cols=2, program_name='Tm calculator(Celsius)', default_size=(700, 740), header_bg_color= '#DCDCDC', terminal_font_color= '#000000', terminal_panel_color= '#DCDCDC',terminal_font_size=14, clear_before_run=True)
def main():
	ap = GooeyParser()
	ap.add_argument("-forward", "--Forward primer", required=True, type=str, help="forward primer sequence")
	ap.add_argument("-reverse", "--Reverse primer", required=True, type=str, help="reverse primer sequence")
	ap.add_argument("-Na", "--Na", required=False, type=float, default=50.0, help="Na+ concentration")
	ap.add_argument("-Mg", "--Mg", required=False, type=float, default=1.5, help="Mg2+ concentration")
	ap.add_argument("-dNTPs", "--dNTPs", required=False, type=float, default=0.6, help="dNTPs concentration")
	ap.add_argument("-additive", "--additive", required=False, type=str, default="no", choices=["no","yes"], help="correct Tm calculation based on additive/s concertrations")
	ap.add_argument("-type", "--additive type", required=False, type=str, default="DMSO", choices=["DMSO","formamide"], help="additive type used")
	ap.add_argument("-per", "--percentage of DSMO", required=False, type=float, default=2.0, help="%DMSO concertration")
	ap.add_argument("-decrease", "--decrease of DMSO", required=False, type=float, default=0.65, help="How much should Tm decrease per percent DMSO")
	ap.add_argument("-conc", "--percentage of formamide", required=False, type=float, default=1.0, help="".join(["%","formamide concertration"]))
	ap.add_argument("-dec", "--decrease of formamide", required=False, type=float, default=0.65, help="How much should Tm decrease per percent formamide")	
	args = vars(ap.parse_args())
	# import forward primer
	frwd = Seq(args['Forward primer'])
	# import reverse primer
	rever = Seq(args['Reverse primer'])
	# calculate Tm of forward primer
	tmf = round(mt.Tm_NN(frwd, Na=args['Na'], Mg=args['Mg'], dNTPs=args['dNTPs']), 2)
	# calculate Tm of reverse primer
	tmr = round(mt.Tm_NN(rever, Na=args['Na'], Mg=args['Mg'], dNTPs=args['dNTPs']), 2)
	# if DMSO or formamide is added fix Tm calculation
	if args['additive'] == "no":
		print(f'The Tm of the forward primer is: {tmf}') # creates the celsius symbol
		print(f'The Tm of the reverse primer is: {tmr}')
		print(f'The Tm difference between the primers is: {abs(round(tmf-tmr, 2))}') # abs will keep the difference positive
	else:
		if args['additive type'] == "DMSO":
			fixedtmf = round(mt.chem_correction(tmf, DMSO=args['percentage of DSMO'],DMSOfactor=args['decrease of DMSO']) ,2)
			print(f'The Tm of the forward primer is: {fixedtmf}')
			fixedtmr = round(mt.chem_correction(tmr, DMSO=args['percentage of DSMO'],DMSOfactor=args['decrease of DMSO']) ,2)
			print(f'The Tm of the reverse primer is: {fixedtmr}')
			print(f'The Tm difference between the primers is: {abs(round(fixedtmf-fixedtmr, 2))}') # abs will keep the difference positive
		else:
			fixedtmf = round(mt.chem_correction(tmf, fmd=args['percentage of formamide'],fmdfactor=args['decrease of formamide']) ,2)
			print(f'The Tm of the forward primer is: {fixedtmf}')
			fixedtmr = round(mt.chem_correction(tmr, fmd=args['percentage of formamide'],fmdfactor=args['decrease of formamide']) ,2)
			print(f'The Tm of the reverse primer is: {fixedtmr}')
			print(f'The Tm difference between the primers is: {abs(round(fixedtmf-fixedtmr, 2))}') # abs will keep the difference positive

if __name__ == '__main__':
	main()

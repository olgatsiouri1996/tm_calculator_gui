# python3
from gooey import *
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
# imput parameters
@Gooey(required_cols=2, program_name='Tm calculator(Celsius)', default_size=(740, 740), header_bg_color= '#DCDCDC', terminal_font_color= '#000000', terminal_panel_color= '#DCDCDC',terminal_font_size=20, clear_before_run=True)
def main():
	ap = GooeyParser()
	ap.add_argument("-forward", "--Forward primer", required=True, type=str, help="Forward primer sequence")
	ap.add_argument("-reverse", "--Reverse primer", required=True, type=str, help="Reverse primer sequence")
	ap.add_argument("-Na", "--Na", required=False, type=float, default=50.0, help="Na+ concentration")
	ap.add_argument("-Mg", "--Mg", required=False, type=float, default=1.5, help="Mg2+ concentration")
	ap.add_argument("-dNTPs", "--dNTPs", required=False, type=float, default=0.6, help="dNTPs concentration")
	ap.add_argument("-additive", "--PCR additive", required=False, type=str, default="no", choices=["no","yes"], help="Correct Tm calculation based on the additive's concertrations")
	ap.add_argument("-type", "--Additive type", required=False, type=str, default="DMSO", choices=["DMSO","formamide"], help="Additive type used")
	ap.add_argument("-per", "--Percentage of DSMO", required=False, type=float, default=2.0, help="%DMSO concertration")
	ap.add_argument("-Decrease", "--Decrease of DMSO", required=False, type=float, default=0.65, help="How much should Tm Decrease per percent DMSO")
	ap.add_argument("-conc", "--Percentage of formamide", required=False, type=float, default=1.0, help="".join(["%","formamide concertration"]))
	ap.add_argument("-dec", "--Decrease of formamide", required=False, type=float, default=0.65, help="How much should Tm Decrease per percent formamide")	
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
	if args['PCR additive'] == "no":
		# print to GUI console
		print(f'The Tm of the forward primer is: {tmf}')
		print(f'The Tm of the reverse primer is: {tmr}')
		print(f'The Tm difference between the primers is: {abs(round(tmf-tmr, 2))}') # abs will keep the difference positive
	else:
		if args['Additive type'] == "DMSO":
			fixedtmf = round(mt.chem_correction(tmf, DMSO=args['Percentage of DSMO'],DMSOfactor=args['Decrease of DMSO']) ,2)
			fixedtmr = round(mt.chem_correction(tmr, DMSO=args['Percentage of DSMO'],DMSOfactor=args['Decrease of DMSO']) ,2)
		else:
			fixedtmf = round(mt.chem_correction(tmf, fmd=args['Percentage of formamide'],fmdfactor=args['Decrease of formamide']) ,2)
			fixedtmr = round(mt.chem_correction(tmr, fmd=args['Percentage of formamide'],fmdfactor=args['Decrease of formamide']) ,2)
		# print to GUI console
		print(f'The Tm of the forward primer is: {fixedtmf}')
		print(f'The Tm of the reverse primer is: {fixedtmr}')
		print(f'The Tm difference between the primers is: {abs(round(fixedtmf-fixedtmr, 2))}') # abs will keep the difference positive

if __name__ == '__main__':
	main()

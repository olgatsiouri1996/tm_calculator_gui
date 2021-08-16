# python3
from gooey import *
# imput parameters
@Gooey(required_cols=2, program_name='Tm calculator', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-forward", "--forward", required=True, type=str, help="forward primer sequence")
	ap.add_argument("-reverse", "--reverse", required=True, type=str, help="forward primer sequence")
	args = vars(ap.parse_args())
# create the function to calculate the Tm:
	def tm(x):
		return 2 * x.count("A") + 2 * x.count("C") + 4 * x.count("G") + 4 * x.count("C")
# main
# import forward primer
	frwd = args['forward']
# import reverse primer
	rever = args['reverse']
# calculate Tm of forward primer
	print("The Tm of the forward primer is:",tm(frwd),sep=" ")
# calculate Tm of reverse primer
	print("The Tm of the reverse primer is:",tm(rever),sep=" ")

if __name__ == '__main__':
	main()

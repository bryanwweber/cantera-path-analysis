import argparse as arg
from .path_analysis import PathAnalysis


def main():
    parser = arg.ArgumentParser(
        description="Create a path analysis file from a CanSen output file"
        )
    parser.add_argument('-c', '--chem', help='Specify the chemistry file')
    parser.add_argument('-s', '--save', help='Specify the stored save file', default='save.hdf')
    parser.add_argument('-o', '--output', help='Specify the output file', default='output.xlsx')
    parser.add_argument('-p', '--percent', help='Specify the integration percent',
                        default=20.0, type=float)
    parser.add_argument('fuel', help='Specify the fuel in the mixture')
    args = parser.parse_args()

    PathAnalysis(
        save=args.save,
        output=args.output,
        percent=args.percent,
        fuel=args.fuel,
        chem=args.chem,
        )

if __name__ == '__main__':
    main()

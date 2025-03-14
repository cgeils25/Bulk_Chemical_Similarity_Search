import argparse

def print_args(args: argparse.Namespace) -> None:
    """
    Print the arguments passed to the script
    """
    print('-'*100)
    print(f'Arguments passed to {__file__}:')
    for arg, value in vars(args).items():
        print(f'{arg}: {value}')
    print('-'*100)

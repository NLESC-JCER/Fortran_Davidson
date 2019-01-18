import argparse

msg = "namd.py -i input"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True, help="program to run")


def main():
    executable = read_cmd_line()
    print("path to executable: ", executable)


def read_cmd_line():
    """
    Read the file to run
    """
    args = parser.parse_args()
    return args.i


if __name__ == "__main__":
    main()

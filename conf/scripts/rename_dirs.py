
def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', type=str, help='Main directory with items that need to be renamed.')
    args = parser.parse_args()

    def rename(entry, replacing=""):

        new_name = entry.name.replace(replacing, "")
        new_path = os.path.join(os.path.dirname(entry.path), new_name)

        print(entry.name, new_name, new_path)

        os.rename(entry.path, new_path)

        return new_path

    for entry in os.scandir(args.dir):
        if entry.name.startswith("job-"):
            rename(entry=entry, replacing="job-")

        elif entry.name.startswith("proj-"):

            new_path = rename(entry=entry, replacing="proj-")

            # Go one level deeper in project directories to rename the data folders.
            for data in os.scandir(new_path):
                rename(entry=data, replacing="store-")


if __name__ == "__main__":
    main()

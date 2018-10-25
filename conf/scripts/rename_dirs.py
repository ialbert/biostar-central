
def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('job_dir', type=str, help='Jobs directory used to rename')
    args = parser.parse_args()

    for entry in os.scandir(args.job_dir):
        if entry.name.startswith("job-"):
            new_name = entry.name.replace("job-", "")
            new_path = os.path.join(os.path.dirname(entry.path), new_name)
            print (entry.name, new_name,new_path)
            os.rename(entry.path, new_path)


if __name__ == "__main__":
    main()

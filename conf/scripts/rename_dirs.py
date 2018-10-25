
def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('job_dir', type=str, help='Jobs directory used to rename')
    args = parser.parse_args()

    for job in os.scandir(args.job_dir):
        if job.name.startswith("job-"):
            new_name = job.name.replace("job-", "")
            new_path = os.path.join(os.path.dirname(job.path), new_name)

            os.rename(job.path, new_path)


if __name__ == "__main__":
    main()

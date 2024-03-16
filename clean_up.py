#!/usr/bin/env python

import subprocess
import shutil
import os
import pprint


def getLast():
    run = subprocess.run(
        "nextflow log | tail -n 1 | cut -f 3",
        shell=True,
        capture_output=True,
    )
    if not run.returncode:
        return run.stdout.decode()


def logRun(run):
    run = subprocess.run(
        f"nextflow log {run}",
        shell=True,
        capture_output=True,
    )
    if not run.returncode:
        return run.stdout.decode().splitlines()


def getOtherDirs():
    run = subprocess.run(
        "find ./work/ -maxdepth 2 -mindepth 2 -type d -exec readlink -f {} \\;",
        shell=True,
        capture_output=True,
    )
    if not run.returncode:
        return run.stdout.decode().splitlines()


def getEmpty():
    run = subprocess.run(
        "find ./work/ -type d -empty -exec readlink -f {} \\;",
        shell=True,
        capture_output=True,
    )
    if not run.returncode:
        return run.stdout.decode().splitlines()


def killWork(work_dirs):
    for f in work_dirs:
        if os.path.isdir(f):
            print(f"Deleting {f}")
            shutil.rmtree(f)


def cli():
    last_run = getLast()
    keep = logRun(last_run)
    others = getOtherDirs()
    while True:
        ask = input(
            "\n"
            "   What would you like to do? [V]iew last run name, [L]ist last "
            "run paths, [C]lean up paths, [E]xit "
        ).lower()
        if ask == "e":
            break
        elif ask == "v":
            print(f"Last run: {last_run}")
        elif ask == "l":
            pprint.pprint(keep)
        elif ask == "c":
            confirm = input("Run nextflow clean? [Y/N] ").lower()
            if confirm == "y":
                subprocess.run(
                    f"nextflow clean -but {last_run.strip()} -f",
                    capture_output=True,
                    shell=True,
                )
            to_delete = set(others) - set(last_run)
            print("The following work directories will be deleted")
            pprint.pprint(to_delete)
            confirm = input("Continue? [Y/N]").lower()
            if confirm == "y":
                killWork(to_delete)
                killWork(getEmpty())


if __name__ == "__main__":
    cli()

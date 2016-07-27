import gzip

# This is used to open files, and if they are .gz files, use gzip to open them.
# Otherwise, open them normally.
def openFile(fileName, *args):
    if fileName[-3:] == ".gz":
        # If there aren't any arguments, default to "rt" - read text. Otherwise 
        # gzip uses "rb" - read binary, and that is not what we want.
        return gzip.open(fileName, *args if len(args) else "rt")
    else:
        return open(fileName, *args)


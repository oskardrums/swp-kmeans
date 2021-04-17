from invoke import task, call

@task
def run(c, n=0, k=0, Random=True):
    c.run("python3 setup.py build_ext --inplace")
    if Random:
        c.run("python3 main.py -r", pty=True)
    else:
        c.run("python3 main.py -n %lu -k %lu" % (n, k), pty=True)


@task(aliases=['del'])
def delete(c):
    c.run("rm -r *.so build/")

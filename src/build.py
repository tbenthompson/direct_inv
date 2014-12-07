from fabricate import *

test_sources = ['test', 'eval']
bench_sources = ['bench', 'eval']

cpp_flags = '-Wall -std=c++11 -Og -DDEBUG'.split()
cpp_flags.append('-I./eigen')
link_flags = []

def build():
    compile()
    link()

def compile():
    for source in test_sources:
        run('g++', '-c', source+'.cpp', cpp_flags)
    for source in bench_sources:
        run('g++', '-c', source+'.cpp', cpp_flags)

def link():
    test_objects = [s+'.o' for s in test_sources]
    bench_objects = [s+'.o' for s in bench_sources]
    run('g++', '-o', 'test', test_objects)
    run('g++', '-o', 'bench', bench_objects)

def clean():
    autoclean()

main()

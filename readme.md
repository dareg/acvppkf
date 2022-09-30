You might need to set the stack to higher value to run the code.
`
ulimit -s unlimited
export OMP_STACKSIZE=4g
`

Enable timers with `export LHOOK=1`. Results can be grouped with the parse_time.awk script.

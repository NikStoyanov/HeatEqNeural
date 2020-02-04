# TODO
- Try to increase ADAM alpha.
- Check if GPU can be used.
- Convergence check on mesh. This is not really a big deal since the process is to be rough
until a closer initial guess is made and then repeat.
- Make figures presentable.
- Save figures and `p` at every 50-100 iterations.
- `p` must be read and saved to a file at every iterations to recover work as a `csv` file
- Exit optimization when $L \leq \sum_{t_{max}^{t_{min}}} \sigma$ where $\sigma$ is TODO.
This is necessary since the optimization will start to capture noise.
- Run on CSF.

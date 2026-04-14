def embed(glob, loc):
    import os
    from ptpython.repl import embed, run_config
    import pathlib

    def configure(repl):
        path = os.path.expanduser('~/.config/ptpython/config.py')
        config_path = pathlib.Path(path)
        run_config(repl, config_path)

    history_path = pathlib.Path(
        os.path.expanduser("~/.local/share/ptpython/history")
    )

    embed(
        glob,
        loc,
        vi_mode=True,
        configure=configure,
        history_filename=history_path,
    )

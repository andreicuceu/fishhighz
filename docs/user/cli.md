# Command line

The installed entry point and Python module call the same `Forecast` interface:

```bash
fishhighz-forecast my-survey.ini --output my-forecast
python -m fishhighz.cli my-survey.ini --output my-forecast
```

Omit the optional positional INI to use the bundled DESI Run-2 recipe:

```bash
fishhighz-forecast --output accuracy-desi2
```

| Argument | Requirement | Meaning |
| --- | --- | --- |
| `ini` | Optional | Native INI path; omitted means the bundled recipe. |
| `--output` | Required | New result directory. Existing directories are rejected. |
| `-h`, `--help` | Optional | Display usage and exit. |

The command prepares and runs a real forecast, then writes the
[JSON and NPZ results](results.md). It has no dry-run, profile, overwrite or
numerical-control flags and no `[output]` INI section. Use the
[INI options](ini.md) for survey and integration settings, or the
[Python interface](python.md) for injection and `run()` arguments. To validate
syntax without a forecast, call `parse_survey_ini()` from Python.

Relative external resource paths are resolved against the INI directory.
`package:` resources resolve inside the installed package. The output path is
relative to the command's current directory. The native calculation needs the
`camb`, `templates` and `survey` extras described in
[Getting started](../getting-started.md).

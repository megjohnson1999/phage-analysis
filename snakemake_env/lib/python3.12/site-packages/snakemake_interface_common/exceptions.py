__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import sys
import textwrap
from typing import Optional

from snakemake_interface_common.rules import RuleInterface


class ApiError(Exception):
    pass


class WorkflowError(Exception):
    def format_arg(self, arg):
        if isinstance(arg, str):
            return arg
        elif isinstance(arg, WorkflowError):
            spec = self._get_spec(arg)

            if spec:
                spec = f" ({spec})"

            return "{}{}:\n{}".format(
                arg.__class__.__name__, spec, textwrap.indent(str(arg), "    ")
            )
        elif sys.version_info >= (3, 11) and isinstance(
            arg,
            ExceptionGroup,  # noqa: F821
        ):
            return "\n".join(self.format_arg(exc) for exc in arg.exceptions)
        else:
            return f"{arg.__class__.__name__}: {arg}"

    def __init__(
        self,
        *args,
        lineno: Optional[int] = None,
        snakefile: Optional[Path] = None,
        rule: Optional[RuleInterface] = None,
    ):
        if rule is not None:
            self.lineno = rule.lineno
            self.snakefile = rule.snakefile
        else:
            self.lineno = lineno
            self.snakefile = snakefile
        self.rule = rule

        # if there is an initial message, append the spec
        if args and isinstance(args[0], str):
            spec = self._get_spec(self)
            if spec:
                args = [f"{args[0]} ({spec})"] + list(args[1:])

        super().__init__("\n".join(self.format_arg(arg) for arg in args))

    @classmethod
    def _get_spec(cls, exc):
        spec = ""
        if exc.rule is not None:
            spec += f"rule {exc.rule.name}"
        if exc.snakefile is not None:
            if spec:
                spec += ", "
            spec += f"line {exc.lineno}, {exc.snakefile}"
        return spec


class InvalidPluginException(ApiError):
    def __init__(self, plugin_name: str, message: str):
        super().__init__(f"Error loading Snakemake plugin {plugin_name}: {message}")

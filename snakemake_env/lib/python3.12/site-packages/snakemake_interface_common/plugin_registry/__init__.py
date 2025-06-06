__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster, Vanessa Sochat"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from abc import ABC, abstractmethod
import types
import pkgutil
import importlib
from typing import List, Mapping

from snakemake_interface_common.exceptions import InvalidPluginException
from snakemake_interface_common.plugin_registry.plugin import PluginBase
from snakemake_interface_common.plugin_registry.attribute_types import AttributeType


class PluginRegistryBase(ABC):
    """This class is a singleton that holds all registered executor plugins."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if hasattr(self, "plugins"):
            # init has been called before
            return
        self.collect_plugins()

    def get_registered_plugins(self) -> List[str]:
        """Return a list of registered plugin names."""
        return [name for name in self.plugins.keys()]

    def is_installed(self, plugin_name: str) -> bool:
        """Return True if the plugin is registered."""
        return plugin_name in self.plugins

    def get_plugin(self, plugin_name: str) -> PluginBase:
        """Get a plugin by name."""
        try:
            return self.plugins[plugin_name]
        except KeyError:
            raise InvalidPluginException(
                plugin_name,
                f"The package {self.module_prefix.replace('_', '-')}{plugin_name} is "
                "not installed.",
            )

    def get_plugin_package_name(self, plugin_name: str) -> str:
        """Get the package name of a plugin by name."""
        return f"{self.module_prefix.replace('_', '-')}{plugin_name}"

    def register_cli_args(self, argparser):
        """Add arguments derived from self.executor_settings to given
        argparser."""
        for _, plugin in self.plugins.items():
            plugin.register_cli_args(argparser)

    def collect_plugins(self):
        """Collect plugins and call register_plugin for each."""
        self.plugins = dict()

        # Executor plugins are externally installed plugins named
        # "snakemake_executor_<name>".
        # They should follow the same convention if on pip,
        # snakemake-executor-<name>.
        # Note that these will not be detected installed in editable
        # mode (pip install -e .).
        for moduleinfo in pkgutil.iter_modules():
            if not moduleinfo.ispkg or not moduleinfo.name.startswith(
                self.module_prefix
            ):
                continue
            module = importlib.import_module(moduleinfo.name)
            self.register_plugin(moduleinfo.name, module)

    def register_plugin(self, name: str, plugin: types.ModuleType):
        """Validate and register a plugin.

        Does nothing if the plugin is already registered.
        """
        if name in self.plugins:
            return

        self.validate_plugin(name, plugin)

        # Derive the shortened name for future access
        plugin_name = name.removeprefix(self.module_prefix).replace("_", "-")

        self.plugins[plugin_name] = self.load_plugin(plugin_name, plugin)

    def validate_plugin(self, name: str, module: types.ModuleType):
        """Validate a plugin for attributes and naming"""
        expected_attributes = self.expected_attributes()
        for attr, attr_type in expected_attributes.items():
            # check if attr is missing and fail if it is not optional
            if not hasattr(module, attr):
                if attr_type.is_optional:
                    continue
                raise InvalidPluginException(name, f"plugin does not define {attr}.")

            attr_value = getattr(module, attr)
            if attr_type.is_class:
                # check for class type
                if not issubclass(attr_value, attr_type.cls):
                    raise InvalidPluginException(
                        name,
                        f"{attr} must be a subclass of "
                        f"{attr_type.cls.__module__}.{attr_type.cls.__name__}.",
                    )
            else:
                # check for instance type
                if not isinstance(attr_value, attr_type.cls):
                    raise InvalidPluginException(
                        name, f"{attr} must be of type {attr_type.cls.__name__}."
                    )

    @property
    @abstractmethod
    def module_prefix(self) -> str: ...

    @abstractmethod
    def load_plugin(self, name: str, module: types.ModuleType) -> PluginBase:
        """Load a plugin by name."""
        ...

    @abstractmethod
    def expected_attributes(self) -> Mapping[str, AttributeType]: ...

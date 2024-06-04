from typing import Dict
from py_script_template.cache import CacheHandler
from py_script_template.cli import get_cli_argument
from py_script_template.config import ConfigHandler, load_config_from_toml


class BaseProvider:
    def __init__(self, config_path: str, cli_config_path: str, cache_path: str):
        self._config_handler: ConfigHandler = load_config_from_toml(
            config_path=config_path
        )
        self._cil_args: Dict = get_cli_argument(cli_config_path)
        self._cache_handler: CacheHandler = CacheHandler(cache_path)

    @property
    def config_handler(self) -> ConfigHandler:
        return self._config_handler

    @config_handler.setter
    def config_handler(self, v: ConfigHandler) -> None:
        self._config_handler = v

    @config_handler.setter
    def config_handler(self, v: ConfigHandler) -> None:
        self._config_handler = v

    @property
    def cache_handler(self) -> CacheHandler:
        return self._cache_handler

    @property
    def cil_args(self) -> Dict:
        return self._cil_args


# class AppProvider(BaseProvider):
#     pass

global_provider = BaseProvider("config.toml", "cli_config.toml", "./cache.json")

if __name__ == "__main__":
    from py_script_template.config import ConfigFields

    handler = global_provider.config_handler
    print(handler.get(ConfigFields.LOG_FILE_SAVE_PATH, "Null Config"))
    print(global_provider.cil_args)

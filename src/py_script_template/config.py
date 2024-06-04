from pydantic import BaseModel
from typing import Any, Dict
import toml


# 这里的这种映射写法其实是为了开发时候直接方便调用已有内容
# 这里的实现类似于ENUM类型 Like This: ConfigFields.NAME
class ConfigFields:
    LOG_LEVEL = "Log_level"
    LOG_CONSOLE_LEVEL = "Log_console_level"
    LOG_FILE_LEVEL = "Log_file_level"
    LOG_FILE_SAVE_PATH = "Log_file_save_path"
    CACHE_FILE_SAVE_PATH = "Cache_file_save_path"
    RESULT_FILE_SAVE_DIR = "Result_file_save_dir"


class ConfigHandler(BaseModel):
    Log_level: int | None
    Log_console_level: int | None
    Log_file_level: int | None
    Log_file_save_path: str | None
    Cache_file_save_path: str | None
    Result_file_save_dir: str | None

    def get(self, name: str, defaultValue: Any = None) -> Any:
        return self.model_dump().get(name, defaultValue)

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    def items(self) -> list[tuple[str, Any]]:
        return list(self.model_dump().items())


def load_config_from_toml(config_path: str) -> ConfigHandler:
    # 读取TOML配置文件
    with open(config_path, "r", encoding="utf-8") as f:
        config_data = toml.load(f)

    # 创建ConfigModel实例并返回
    return ConfigHandler(**config_data)


# 调用示例
if __name__ == "__main__":
    config_file = "./config.toml"  # 替换为实际的配置文件路径
    config = load_config_from_toml(config_file)
    print(config.model_dump())

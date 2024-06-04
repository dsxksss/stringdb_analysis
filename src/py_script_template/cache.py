# TODO 待使用tinydb重构

from typing import Any, Dict, List
from tinydb import TinyDB, Query


class CacheHandler:
    def __init__(self, save_path: str):
        self.db = TinyDB(save_path)

    # 检查数据是否已存在
    def is_exists(self, name: str) -> bool:
        result = self.db.search(Query().name == name)
        return len(result) > 0

    # 存储缓存
    def save_cache(self, name: str, key: str, value: Any) -> None:
        if not self.is_exists(name):
            self.db.insert({"name": name, key: value})

    # 更新缓存
    def update_cache(self, name: str, key: str, value: Any) -> None:
        self.db.update({key: value}, Query().name == name)

    # 获取缓存
    def get_cache(self, name: str) -> Dict[str, Any]:
        data = self.db.search(Query().name == name)
        result = data[0]
        return result if result else dict()

    # 获取全部缓存
    def get_cache_list(self) -> List[Any]:
        return self.db.all()

    # 根据缓存名删除某个缓存
    def remove_cache(self, name: str) -> None:
        self.db.remove(Query().name == name)

    # 清空缓存
    def clear_all_cache(self) -> None:
        self.db.clear_cache()

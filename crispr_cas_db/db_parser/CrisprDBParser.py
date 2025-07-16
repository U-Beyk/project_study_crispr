'''
This parser will load the sql dump of the CRISPR CAS Database and will process it. 
It will create a json-file for every table inside the sql dump.

author: U.B.
'''

from re import Match
import re
import json
import os

class CrisprDBParser:
    '''Class for parsing the sql dump of the CRISPR CAS database'''
    def __init__(self, sql_file_path: str = "./crispr_cas_db/db_parser/Crispr_Cas_Database_SQL_Dump.sql") -> None:
        '''Initializes the CrisprDBParser object'''
        self._sql_file_path = sql_file_path
        self._sql_table_pattern = re.compile(r'^COPY\s+public.(\S+)\s*\((.*?)\)')
        self._in_table = False
        self._table = None

    def process_sql_file(self) -> None:
        '''Processes the sql dump line by line and creates json files for eavh table'''
        with open(self._sql_file_path, "r") as file:
            for line in file:
                self._process_line(line)

    def _process_line(self, line: str) -> None:
        '''Processes each line by checking for certain regular expressions and handles them accordingly'''
        match = self._sql_table_pattern.match(line)
        if match:
            self._table = CrisprTable(match)
            self._in_table = True
            return
        if self._in_table:
            self._process_sql_table(line.strip())
    
    def _process_sql_table(self, stripped_line: str) -> None:
        '''Processes the sql table by either adding an entry or saving the finished table as a json'''
        if not stripped_line:
            self._in_table = False
            self._save_json()
        else:
            self._table.add_entry(stripped_line)

    def _save_json(self) -> None:
        '''Saves the finished table as a json'''
        os.makedirs("crispr_cas_db/database_tables", exist_ok=True)
        json_file_name = f"{self._table.table_name}.json"
        with open(f"crispr_cas_db/database_tables/{json_file_name}", "w") as json_file:
            json.dump(self._table.table_content, json_file, indent=4)
        print(f"Json-file: {json_file_name} was created")


class CrisprTable:
    '''A class representing a table from the CRISPR CAS database'''
    def __init__(self, pattern_match: Match) -> None:
        '''Initializes a CrisprTable object'''
        self._table_name = pattern_match.group(1)
        self._table_columns = [f"{self._table_name}_{column}" for column in pattern_match.group(2).split(", ")]
        self._table_content = []

    @property
    def table_name(self) -> str:
        '''Gets and returns the table name as a string'''
        return self._table_name
    
    @property
    def table_content(self) -> list[dict[str, str]]:
        '''Gets and returns the table content as a list'''
        return self._table_content

    def add_entry(self, entry: str) -> None:
        '''Adds an entry from the sql dump'''
        entry_list = entry.split("\t")
        entry_dictionary = dict(zip(self._table_columns, entry_list))
        self._table_content.append(entry_dictionary)
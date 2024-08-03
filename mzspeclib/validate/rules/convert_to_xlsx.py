import argparse
import itertools
import json

import csv
import xlsxwriter

import psims

cv = psims.load_psims()

def rule_to_table_rows(rule: dict):
    row_base = {
        "id": rule["id"],
        "path": rule["path"],
        "level": rule.get("level", rule.get("requirement_level")),
        "notes": rule.get("notes"),
        "combination_logic": rule.get("combination_logic"),
    }

    if "condition" in rule:
        condition = rule["condition"]
        row_base["condition_name"] = condition["name"]
        row_base["condition_accession"] = condition["accession"]
        row_base["condition_allow_children"] = condition["allow_children"]
        if "value" in condition:
            cond_value = condition["value"]
            for k, v in cond_value.items():
                row_base[f"condition_value_{k}"] = k
                row_base[f"condition_value_{v}"] = v

    if "attr" not in rule:
        yield row_base
    else:
        for attr in rule["attr"]:
            row = row_base.copy()
            row["attribute_name"] = attr["name"]
            row["attribute_accession"] = attr["accession"]
            row["allow_children"] = attr.get("allow_children")
            row['repeatable'] = attr.get('repeatable', False)
            term = cv[attr["accession"]]
            row["description"] = term.definition
            units = ";".join([f"{rel.accession}|{rel.comment}" for rel in term.get("has_units", [])])
            if units:
                row["has_units"] = units

            row['default_unit'] = attr.get("default_unit")
            row["attr_notes"] = attr.get("notes")
            if "value" in attr:
                value_rule = attr["value"]
                if isinstance(value_rule, str):
                    row["value_rule"] = value_rule
                else:
                    row["value_rule"] = value_rule["name"]
                    value_constraint = None
                    for k, v in value_rule.items():
                        if k != "name":
                            if value_constraint is not None:
                                raise ValueError("Too many properties in value rule")
                            else:
                                value_constraint = v
                    row["value_constraint"] = value_constraint
            if "condition" in attr:
                condition = attr["condition"]
                row["condition_name"] = condition["name"]
                row["condition_accession"] = condition["accession"]
                row["condition_allow_children"] = condition["allow_children"]
                if "value" in condition:
                    cond_value = condition["value"]
                    for k, v in cond_value.items():
                        row[f"condition_value_{k}"] = k
                        row[f"condition_value_{v}"] = v

            yield row


def infer_headers(rows):
    seen = set(("notes",))
    keys = []
    for row in rows:
        for k in row:
            if k in seen:
                continue
            seen.add(k)
            keys.append(k)
    keys.append("notes")
    return keys


def create_issue_link(row):
    link = f"https://github.com/HUPO-PSI/mzSpecLib/issues/new?labels=validation&title={row['id'].replace(' ', '+')}"
    return link


def rule_set_to_worksheet(rule_set, sheet, wrap_fmt, default_fmt):
    rows = list(itertools.chain.from_iterable(map(rule_to_table_rows, rule_set["rules"])))
    headers = infer_headers(rows)
    sheet.write_row(0, 0, headers + ["link"], default_fmt)

    for i, row in enumerate(rows, 1):
        sheet.write_row(i, 0, [row.get(k) for k in headers])
        sheet.write_url(i, len(headers), create_issue_link(row), string="Open Issue")

    for i, h in enumerate(headers):
        if h == "notes":
            sheet.set_column(i, i, 60, cell_format=wrap_fmt)
        elif h == "attribute_name":
            sheet.set_column(i, i, 30, cell_format=default_fmt)
        elif h == "description":
            sheet.set_column(i, i, 60, cell_format=wrap_fmt)
        elif h in ("level", "combination_logic"):
            sheet.set_column(i, i, 10, cell_format=default_fmt)
        else:
            sheet.set_column(i, i, 20, cell_format=default_fmt)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outpath", help="The name of the file to write")
    parser.add_argument("rule_sets", nargs="+", help="The rule set JSON file paths to convert")
    parser.add_argument("-f", "--format", choices=["xlsx", "csv"])
    args = parser.parse_args()

    rule_set_paths = list(filter(lambda x: not x.endswith("rules-schema.json"), args.rule_sets))

    if args.format == 'xlsx':
        with xlsxwriter.Workbook(args.outpath) as wb:
            wrap_fmt = wb.add_format({"text_wrap": True, "align": "vcenter"})
            default_fmt = wb.add_format({"align": "vcenter"})

            for path in rule_set_paths:
                rule_set = json.load(open(path))
                name = rule_set["name"]
                sheet = wb.add_worksheet(name)
                rule_set_to_worksheet(rule_set, sheet, wrap_fmt, default_fmt)

    elif args.format == 'csv':
        with open(args.outpath, 'wt', newline='') as fh:
            all_rows = []
            for path in rule_set_paths:
                rule_set = json.load(open(path))
                rows = list(itertools.chain.from_iterable(map(rule_to_table_rows, rule_set["rules"])))
                for row in rows:
                    row['rule_set_name'] = rule_set['name']
                    row["link"] = create_issue_link(row)
                all_rows.extend(rows)

            headers = infer_headers(all_rows)
            headers.remove("rule_set_name")
            headers.remove("link")
            headers = ['rule_set_name'] + headers + ['link']
            writer = csv.DictWriter(fh, headers)
            writer.writeheader()
            writer.writerows(all_rows)
            fh.flush()



if __name__ == "__main__":
    main()
#legacy
def add_rs(conn, task):
    sql = '''INSERT INTO dbsnp(id, chromosome, coordinate, rs, refAllele, altAllele, refAlleleFrequency, altAlleleFrequencySum)
             VALUES(?,?,?,?,?,?,?,?) '''
    cur = conn.cursor()
    cur.execute(sql, task)
    conn.commit()
    return cur.lastrowid
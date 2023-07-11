import sqlite3

db = sqlite3.connect("theory.db")
c = db.cursor()

res = c.execute("SELECT * FROM TheoryIndex where ID=400").fetchone()
# Now change the ID
new_res = list(res)
new_res[0] = 399
new_res[36] = "Copy of ID=400 with low precision grids, cfactors from Th 200"

# And insert it
qmarks = ",".join(["?"]*len(new_res))
insert_query = f"insert into TheoryIndex values ({qmarks})"
c.execute(insert_query, new_res)
db.commit()

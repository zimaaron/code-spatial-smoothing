print(toc())

(glue("no print time: {toc()}"))

print(glue("with print time: {toc()[3]}"))

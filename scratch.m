conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')

docs = find(conn, 'states')

val = jsondecode(docs(10).game_states)

val

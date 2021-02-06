#!/usr/bin/env python3

import cgi

# Useful debugging output
import cgitb

cgitb.enable()  # Send errors to browser
# cgitb.enable(display=0, logdir="/path/to/logdir") # log errors to a file
 

# Print the HTML MIME-TYPE header
print ("Content-Type: text/html\n")

# Now grab the content of the form
form = cgi.FieldStorage()


import cgi, os
import cgitb; cgitb.enable()

form = cgi.FieldStorage()

# Get filename here.
fileitem = form['filename']
fileID=form["ID"]

# Test if the file was uploaded
if fileitem.filename:
	# strip leading path from file name to avoid 
	# directory traversal attacks
	#fn = os.path.basename(fileitem.filename)    #if windows
	fn = os.path.basename(fileitem.filename.replace("\\", "/" ))
	if fn.endswith('.pdb'):
		save_user_pdb=open(os.path.join('/serv/www/html_lilian/uploads/' + fn), 'wb')
		save_user_pdb.write(fileitem.file.read())
	else:
		print ("""
			<html>
			<body>
			   <p>%s</p>
			</body>
			</html>
			""") % ("error: file is not a PDB file")

   
print ("""
<html>
<body>
   <p>%s</p>
</body>
</html>
""") % ()

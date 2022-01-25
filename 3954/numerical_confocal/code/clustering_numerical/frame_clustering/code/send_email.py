# Import Python Packages
import smtplib
# Set Global Variables
gmail_user = 'martinaoliver1'
gmail_password = 'martiolisardi'
# Create Email 
mail_from = gmail_user
mail_to = 'mo2016@ic.ac.uk'
mail_subject = 'Hello'
mail_message_body = 'Hello World!'

mail_message = '''\
From: %s
To: %s
Subject: %s
%s
''' % (mail_from, mail_to, mail_subject, mail_message_body)
# Sent Email
server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
server.login(gmail_user, gmail_password)
server.sendmail(mail_from, mail_to, mail_message)
server.close()

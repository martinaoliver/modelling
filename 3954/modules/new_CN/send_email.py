import sys
import smtplib
def sendemail():
    gmail_user = 'martioli333@gmail.com'
    gmail_password = 'martiolisardi'
    scriptname = sys.argv[0]
    sent_from = gmail_user
    to = ['mo2016@ic.ac.uk']
    subject = 'Job finished: %s'%str(scriptname)
    body = 'Job finished: %s'%str(scriptname)

    email_text = """\
    From: %s
    To: %s
    Subject: %s

    %s
    """ % (sent_from, ", ".join(to), subject, body)

    try:
        smtp_server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
        smtp_server.ehlo()
        smtp_server.login(gmail_user, gmail_password)
        smtp_server.sendmail(sent_from, to, email_text)
        smtp_server.close()
        print ("Email sent successfully!")
    except Exception as ex:
        print ("Something went wrong….",ex)

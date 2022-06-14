import smtplib  
def sendmail(jobID):
    conn =smtplib.SMTP('smtp.mail.yahoo.com',587)  
    type(conn)  
    conn.ehlo()  
    conn.starttls()  
    conn.login('martinaoliver1@yahoo.com','rrttnutjxbfbhgsf')  
    conn.sendmail('martinaoliver1@yahoo.com','mo2016@ic.ac.uk','Subject:%s'%jobID)  
    conn.quit() 
    print('emailsent')


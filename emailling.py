import smtplib
import mimetypes
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.MIMEImage import MIMEImage
from email.MIMEAudio import MIMEAudio
from email.Utils import COMMASPACE, formatdate
from email import Encoders

from email_text_ITA import welcome_word, email_subject, email_body

def emailling(user_name, From, To, PWD , FilePath, FileNames):
    
    msg = MIMEMultipart()
    msg['From'] = From
    msg['To'] = To
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = email_subject

    msg.attach(MIMEText(welcome_word + user_name))
    msg.attach(MIMEText(email_body))

    try:
        smtp = smtplib.SMTP('smtp.gmail.com:587')
        smtp.starttls()
        smtp.login(From, PWD)
    except:
        i = 1
    else:
        i = 0

    if i == 0:
        for FileName in FileNames:
            file = FilePath + "/" + FileName
            ctype, encoding = mimetypes.guess_type(file)
            if ctype is None or encoding is not None:
                # No guess could be made, or the file is encoded (compressed), so
                # use a generic bag-of-bits type.
                ctype = 'application/octet-stream'
            maintype, subtype = ctype.split('/', 1)
            if maintype == 'text':
                fp = open(file)
                # Note: we should handle calculating the charset
                part = MIMEText(fp.read(), _subtype=subtype)
                fp.close()
            elif maintype == 'image':
                fp = open(file, 'rb')
                part = MIMEImage(fp.read(), _subtype=subtype)
                fp.close()
            elif maintype == 'audio':
                fp = open(file, 'rb')
                part = MIMEAudio(fp.read(), _subtype=subtype)
                fp.close()
            else:
                fp = open(file, 'rb')
                part = MIMEBase(maintype, subtype)
                part.set_payload(fp.read())
                fp.close()
                # Encode the payload using Base64
                Encoders.encode_base64(part)
            part.add_header('Content-Disposition', 'attachment; filename="%s"' %FileName)
            msg.attach(part)
        try:
            smtp.sendmail(From, To, msg.as_string())
        except:
            print "Mail not sent"
        else:
            print "Mail sent"
        smtp.close()
    else:
        print "Connection failed"


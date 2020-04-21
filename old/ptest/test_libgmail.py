import libgmail

ga = libgmail.GmailAccount("erin.sheldon@gmail.com", "tL3NZcUmP")
ga.login()

gc = ga.getContacts()

for contact in gc.contactList:
    print contact.name, "<"+contact.email+">"

#folder = ga.getMessagesByFolder('inbox')

#for thread in folder:
#  print len(thread), thread.subject, thread.snippet
#  #print thread.id, len(thread), thread.subject
#  for msg in thread:
#    print "  ", msg.number, msg.subject
#    #print "  ", msg.id, msg.number, msg.subject
#    #print msg.source

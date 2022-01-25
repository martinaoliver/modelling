#command; echo "Process done" | mail -s "Process done" mo2016@ic.ac.uk
ssmtp mo2016@ic.ac.uk <<<$'Subject: testing 1...2...3'

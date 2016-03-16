#ifndef NETIO_H_
#define NETIO_H_

#include "Bytes.h"
#include <time.h>

class Socket
{
protected:
	int m_socket;

public:
	Socket();
	Socket(int socket) : m_socket(socket) {}
	virtual ~Socket();

	void write_bytes(const Bytes &bytes);
	Bytes read_bytes();

        void write_2_ciphertexts(const Bytes &bytes);
        void write_3_ciphertexts(const Bytes &bytes);
        void write_4_ciphertexts(const Bytes &bytes);
        void write_n_ciphertexts(const Bytes &bytes, const uint32_t n);
        Bytes read_2_ciphertexts();
        Bytes read_3_ciphertexts();
        Bytes read_4_ciphertexts();
        Bytes read_n_ciphertexts(Bytes & bytes, const uint32_t n);

	void write_string(const std::string &str);
	std::string read_string();
};

class ClientSocket : public Socket
{
public:
	ClientSocket(const char *host_ip, size_t port);
	virtual ~ClientSocket() {}
};

class ServerSocket : public Socket
{
	std::vector<int> m_sockets;

public:
	ServerSocket(size_t port);
	Socket *accept();
	virtual ~ServerSocket();
};

inline void my_read(int socket, void *data, size_t n);


#endif /* NETIO_H_ */

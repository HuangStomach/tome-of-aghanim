
package com.example.multidatasource.config;

import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.autoconfigure.orm.jpa.JpaProperties;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.jdbc.DataSourceBuilder;
import org.springframework.boot.orm.jpa.EntityManagerFactoryBuilder;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Primary;
import org.springframework.data.jpa.repository.config.EnableJpaRepositories;
import org.springframework.orm.jpa.JpaTransactionManager;
import org.springframework.orm.jpa.LocalContainerEntityManagerFactoryBean;
import org.springframework.transaction.PlatformTransactionManager;
import org.springframework.transaction.annotation.EnableTransactionManagement;

import javax.persistence.EntityManagerFactory;
import javax.sql.DataSource;
import java.util.HashMap;
import java.util.Map;

@Configuration
@EnableJpaRepositories(
    basePackages = "com.example.multidatasource.repository.b",
    entityManagerFactoryRef = "entityManagerFactoryB",
    transactionManagerRef = "transactionManagerB"
)
public class DataSourceConfigB {
    @Bean(name = "dataSourceB")
    @ConfigurationProperties(prefix = "spring.datasource.b")
    public DataSource dataSourceB() {
        return DataSourceBuilder.create().build();
    }

    @Bean(name = "jpaPropertiesB")
    @ConfigurationProperties(prefix = "spring.jpa.b")
    public JpaProperties jpaPropertiesB() {
        return new JpaProperties();
    }

    @Bean(name = "entityManagerFactoryB")
    public LocalContainerEntityManagerFactoryBean entityManagerFactoryB(
            EntityManagerFactoryBuilder builder,
            @Qualifier("dataSourceB") DataSource dataSource) {
        Map<String, String> props = jpaPropertiesB().getProperties();

        return builder
                .dataSource(dataSource)
                .packages("com.example.multidatasource.entity.b")
                .persistenceUnit("b")
                .properties(props)
                .build();
    }

    @Bean(name = "transactionManagerB")
    public PlatformTransactionManager transactionManagerB(
            @Qualifier("entityManagerFactoryB") EntityManagerFactory entityManagerFactory) {
        return new JpaTransactionManager(entityManagerFactory);
    }
}
